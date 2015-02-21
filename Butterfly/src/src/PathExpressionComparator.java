import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.HashSet;
import java.util.ArrayList;



public class PathExpressionComparator implements Comparator<Object> {


	HashMap<PairPath,HashSet<List<Integer>>> pp_to_transcript;
	HashMap<List<Integer>,Float> transcript_to_fractional_expr;
	HashMap<List<Integer>,Float> transcript_to_fractional_mass;
	HashMap<List<Integer>,HashMap<PairPath,Float>> transcript_to_pp_frag_counts; // split out by pairpath
	HashMap<List<Integer>, Float> transcript_to_sum_frag_counts; // sum of all fractional frag counts
	
	HashMap<PairPath,Integer> pp_to_read_support;
	
	HashMap<List<Integer>, Integer> seqLengths;
	
	public static int MAX_ROUNDS = 100;
	public static int MIN_ROUNDS = 100;
	public static float MIN_DELTA_REQUIRED = 0.01f;
	
	
	
	
	
	public PathExpressionComparator (List<List<Integer>> all_paths, 
									HashMap<List<Integer>,HashMap<PairPath, Integer>> PathReads, 
									HashMap<List<Integer>, Integer> seqLengths) {
		
		pp_to_transcript = new HashMap<PairPath,HashSet<List<Integer>>>();
		transcript_to_fractional_expr = new HashMap<List<Integer>,Float>();
		transcript_to_fractional_mass = new HashMap<List<Integer>,Float>();
		transcript_to_pp_frag_counts = new HashMap<List<Integer>,HashMap<PairPath,Float>>();
		pp_to_read_support = new HashMap<PairPath,Integer>();
		transcript_to_sum_frag_counts = new HashMap<List<Integer>, Float>();
		
		this.seqLengths = seqLengths;
		
		BFLY_GLOBALS.debugMes("EM on: " + all_paths.size() + " paths.", 10);
		
		//if (VERBOSE)
		//	System.err.println("-reorganizing data structures for EM computation.");
		
		// initialize
		for (List<Integer> path : PathReads.keySet()) {
			
			
			// only analyze the paths of interest, not all stored in the more comprehensive PathReads data structure.
			if (! all_paths.contains(path))
				continue;
			
			
			transcript_to_pp_frag_counts.put(path, new HashMap<PairPath,Float>());
			transcript_to_sum_frag_counts.put(path, 0f);
			
			
			HashMap<PairPath,Integer> pp_map = PathReads.get(path);
			pp_to_read_support.putAll(pp_map);
			
			for (PairPath pp : pp_map.keySet()) {
				
				//if (! pp.isCompatibleAndContainedBySinglePath(path))
				//	continue;
				
				if (! pp_to_transcript.containsKey(pp)) {
					pp_to_transcript.put(pp, new HashSet<List<Integer>>());
				}
				pp_to_transcript.get(pp).add(path);
				transcript_to_pp_frag_counts.get(path).put(pp, 0f); // just init for now, store fractional assignments later on.
			}
			
			
		}
		
	
		
		run_EM();
		
		
		
	}
	
	
	
	
	private void run_EM () {
		
		init_transcript_expr();
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
			System.err.println("\n\n============================================\nEM-round[0]");
			describe_expr_values();
		}
		
		boolean enough_iterations = false;
		float prev_likelihood = calc_likelihood();
		
		int round = 0;
		
		while (! enough_iterations) {
			
		
			round++;
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
				System.err.println("\n\n============================================\nEM-round[" + round + "]");
			}
			
		
			
			long start_time = System.currentTimeMillis();
			
			long start_time_E_step = System.currentTimeMillis();
			calc_expected_read_alignments ();
			long end_time_E_step = System.currentTimeMillis();
			//System.err.println("Time for E-step: " + (end_time_E_step - start_time_E_step));
			
			long start_time_M_step = System.currentTimeMillis();
			recompute_expression_values();
			long end_time_M_step = System.currentTimeMillis();
			//System.err.println("Time for M-step: " + (end_time_M_step - start_time_M_step));
			
		
			
			long start_time_likelihood = System.currentTimeMillis();
			float likelihood = calc_likelihood();
			long end_time_likelihood = System.currentTimeMillis();
			//System.err.println("Time for likelihood calc: " + (end_time_likelihood - start_time_likelihood));
			
			float delta = likelihood - prev_likelihood;
			
			long end_time = System.currentTimeMillis();
			
			long seconds_for_computation = (end_time - start_time) / 1000l;
			
			BFLY_GLOBALS.debugMes("EM round[" + round + "] Prev: " + prev_likelihood + ", curr: " + likelihood 
						+  ", delta: " + delta + " [seconds: " + seconds_for_computation + "]\n", 10);
		
			if (delta < -1) {
				//throw new RuntimeException("Error, shouldn't have negative delta during EM: " + delta + ", round: " + round);
				System.err.println("warning, shouldn't have negative delta during EM: " + delta + ", round: " + round);
			}
			
			prev_likelihood = likelihood;
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20)
				describe_expr_values();
			
			
			if ( Math.abs(delta) < MIN_DELTA_REQUIRED || round >= MAX_ROUNDS)
				break; 
		}
		
	
		
	}
	
	
	private void init_transcript_expr() {
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20)
			System.err.println("// initializing EM (round 0)");
		
		init_transcript_to_sum_frags();
		
		
		// loop through pair paths, assign  read_support/num_transcript to each transcript
		for (PairPath pp : pp_to_transcript.keySet()) {
			
			HashSet<List<Integer>> transcript_list = pp_to_transcript.get(pp);
			int num_transcripts_with_pp = transcript_list.size();
			int read_support = pp_to_read_support.get(pp);
			
			// initially, divide multiply mappers among the transcripts
			// including the number of reads in that pp
			float fractional_support = (float)read_support/num_transcripts_with_pp;
			
			for (List<Integer> transcript : transcript_list) {
				
				// add fractional support to the expression value of that transcript.
				
				transcript_to_sum_frag_counts.put(transcript, transcript_to_sum_frag_counts.get(transcript) + fractional_support);
				
				
				// store fractional support for this transcript here.
				this.transcript_to_pp_frag_counts.get(transcript).put(pp, fractional_support);  
				
			}
		}
		
		
		// set initial expression values.
		float sum_expr_values = 0;
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			float frags = transcript_to_sum_frag_counts.get(transcript);
			
			float expr = compute_expr(transcript, frags); 
			
			this.set_expr(transcript, expr); // temporarily set it to frags/length
			
			sum_expr_values += expr;
		}
		
		float sum_mass = 0;
		// convert to relative expression value:
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			float rel_expr_val = get_expr(transcript) / sum_expr_values;
			set_expr(transcript, rel_expr_val); // now set it to fraction all transcripts
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20)
				System.err.println(this.describe_frag_assignment(transcript));
		
			
			float mass = rel_expr_val * seqLengths.get(transcript);
			this.transcript_to_fractional_mass.put(transcript, mass); // temporary, set to fractional value below
			
			sum_mass += mass;
		}
		
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			
			float rel_expr_val = get_expr(transcript);
			
			this.transcript_to_fractional_mass.put(transcript, this.transcript_to_fractional_mass.get(transcript)/sum_mass); // now set to correct value.
			
		}
		
		
	}
	
	
	private void init_transcript_to_sum_frags() {
		
		// initialize all to zero frags
		
		for (List<Integer> path : transcript_to_sum_frag_counts.keySet()) {
			transcript_to_sum_frag_counts.put(path, 0f);
		}
		
	}




	public String describe_frag_assignment (List<Integer> transcript) {
		
		String text = "";
		
		float sum_frags_assigned = 0;
		for (PairPath pp : transcript_to_pp_frag_counts.get(transcript).keySet()) {
			
			float orig_pp_read_support = this.pp_to_read_support.get(pp);

			float frags_assigned = transcript_to_pp_frag_counts.get(transcript).get(pp);
			
			float pct_pp_reads_assigned = frags_assigned / orig_pp_read_support * 100;
			
			text += "FRAGS_ASSIGNED: [" + frags_assigned + ", " + pct_pp_reads_assigned + "%] for read: " + pp + "\n";
			sum_frags_assigned += frags_assigned;
		}
		
		text = "Frag assignment for trans: " + transcript + "\n"
				 + "len=" + this.seqLengths.get(transcript) + " expr=" + get_expr(transcript) + " sum_frags=" + sum_frags_assigned + "\n" + text;
		
		return(text);
		
	}
	
	
	public float get_expr (List<Integer> path) {
		
		if (transcript_to_fractional_expr.containsKey(path)) {
		
			return(transcript_to_fractional_expr.get(path));
		}
		else {
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 10) {
				System.err.println("WARNING: no expr value stored for path: " + path);
			}
			return(0f);
		}
	}
	
	
	private void describe_expr_values () {
		
		System.err.println("####  ACCOUNTING FOR EXPR VALUES");
		
		float sum_expr = 0;
		
		for (List<Integer> transcript : transcript_to_fractional_expr.keySet()) {
			
			float expr = get_expr(transcript);
			sum_expr += expr;
			
			System.err.println("expr: " + expr + " for path: " + transcript);
			
		}
		
		System.err.println("\n** Sum expr: " + sum_expr + "\n\n");	
		
	}
	
	
	
	private Float calc_likelihood () {
		
		//if (VERBOSE)
		//	System.err.println("-computing likelihood");
		
		// compute probability of all observed read alignments, given isoform relative expression values.
		
		//  formula: prod_for_all_reads (sum_all_isoform_mappings_by_read ( isoform_relative_expression *  p(aligned_read) ) )
		//           where p(aligned_read) = (theta_i/length_i)/sum_all(theta_j/length_j)) for those transcripts having the read mapped.
		
		float log_likelihood = 0;
		
		for (PairPath pp : pp_to_transcript.keySet()) {
			HashSet<List<Integer>> transcripts = pp_to_transcript.get(pp);
			
			int count_reads_in_pp = pp_to_read_support.get(pp);
			
			
			float sum_pp_likelihood = 0;
			
			for (List<Integer> transcript : transcripts) {
				//int num_pairpaths_in_transcript = transcript_to_pp_frag_counts.get(transcript).size();
				
				float likelihood = this.transcript_to_fractional_mass.get(transcript) * (1f/this.seqLengths.get(transcript));
				
				sum_pp_likelihood += likelihood;
			
			}
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20)
				System.err.println("sum_pp_log_likelihood: " + (count_reads_in_pp * Math.log(sum_pp_likelihood)) + " for PP: " + pp);
			
			log_likelihood += count_reads_in_pp *  Math.log(sum_pp_likelihood);
			
		}
	
		
		return(log_likelihood);
		
	}
	
	
	//--------------
	// The E-step
	//--------------
	
	private void calc_expected_read_alignments () {
		
		//if (VERBOSE)
		//	System.err.println("-calc expected read alignments (E-step)");
		
		// fractional read assignment
		
		//  Expected value for read alignment =   prob of read alignment to transcript_i /  sum of prob of read aligned to all other transcripts
		
		//                                    =   (expression(transcript_i) / length(transcript_i) ) / sum_for_all_multiply_mapped_read (expr(trans_o)/length(trans_o)))
		
		
		
		this.init_transcript_to_sum_frags();
		
		float sum_all_frags = 0;
		
		//float pp_counter = 0;
		for (PairPath pp : pp_to_transcript.keySet()) {
			HashSet<List<Integer>> transcripts = pp_to_transcript.get(pp);
			
			//pp_counter++;
			//System.err.println("\r[" + pp_counter + "]  ");
			
			
			int count_reads_in_pp = pp_to_read_support.get(pp);
			
			float sum_read_probs = 0;
			for (List<Integer> transcript : transcripts) {
			
				
				float prob = this.transcript_to_fractional_mass.get(transcript) *  (1f/ this.seqLengths.get(transcript));
				
				sum_read_probs += prob; 
				
			}
			
			for (List<Integer> transcript : transcripts) {
				

				float this_read_prob = this.transcript_to_fractional_mass.get(transcript) *  (1f/ this.seqLengths.get(transcript));
				
				float expected_fractional_assignment = this_read_prob / sum_read_probs;
				
				float fractional_frags_assigned = count_reads_in_pp * expected_fractional_assignment;
				
			
				transcript_to_sum_frag_counts.put(transcript, transcript_to_sum_frag_counts.get(transcript) + fractional_frags_assigned);
				
				// store fractional support for this transcript here.
				this.transcript_to_pp_frag_counts.get(transcript).put(pp, fractional_frags_assigned); 
			
			
				sum_all_frags += fractional_frags_assigned;
			}
			
		}
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20)
			System.err.println("Total frags: " + sum_all_frags + "\n"); // verified total mass held constant
		
	}
	
	//------------
	// The M-step
	//------------
	
	private void recompute_expression_values () {
	
		//if (VERBOSE)
		//	System.err.println("-recomputing fractional expression values (M-step)");
		
		
		float sum_expr_vals = 0;
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20)
				System.err.println(describe_frag_assignment(transcript));
			
			float frag_counts = transcript_to_sum_frag_counts.get(transcript);
			//int num_pairpaths_in_transcript = transcript_to_pp.get(transcript).size();
			
			float expr = compute_expr(transcript, frag_counts); // / length of transcript
		
			// temporary replace w/ expr value
			set_expr(transcript, expr);
			
			sum_expr_vals += expr;
		}
		
	
		
		
		// update w/ fractional expression value
		float sum_fractional_expr_vals = 0;
		float sum_mass = 0;
		
		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			
			float expr = get_expr(transcript);
			
			float fractional_expression = expr/sum_expr_vals;
			
			set_expr(transcript, fractional_expression); // now set to the fractional expression value.
				
			sum_fractional_expr_vals += fractional_expression;
		
			float mass = fractional_expression * seqLengths.get(transcript);
			this.transcript_to_fractional_mass.put(transcript, mass); // temporary, set to fractional value below
			
			sum_mass += mass;
		}

		for (List<Integer> transcript : transcript_to_sum_frag_counts.keySet()) {
			

			this.transcript_to_fractional_mass.put(transcript, this.transcript_to_fractional_mass.get(transcript)/sum_mass); // now set to correct value.

		}
		
		//System.err.println("Sum expr vals: " + sum_fractional_expr_vals); // verified should be 1
		
	}
	

	private float compute_expr(List<Integer> transcript, float frag_counts) {
		
		float expr = frag_counts / this.seqLengths.get(transcript);
		
		return(expr);
		
	}
	
	
	private void set_expr(List<Integer> transcript, float fractional_expression) {
		
		this.transcript_to_fractional_expr.put(transcript, fractional_expression);
		
	}




	@Override
	public int compare(Object o1, Object o2) {
		
		List<Integer> path1 = (List<Integer>) o1;
		List<Integer> path2 = (List<Integer>) o2;
		
		
		// Float expr_path1 = this.get_transcript_to_sum_frag_counts(path1); 
		// Float expr_path2 = this.get_transcript_to_sum_frag_counts(path2); 
		
		Float expr_path1 = transcript_to_fractional_expr.get(path1);
		Float expr_path2 = transcript_to_fractional_expr.get(path2);
		
		if (expr_path1 == null)  // not sure why this happens... extremely rare, and on linux not mac
			expr_path1 = 0f;
		
		if (expr_path2 == null)
			expr_path2 = 0f;
		
		
		if (expr_path1 < expr_path2) {
			return(-1);
		}
		else if (expr_path1 > expr_path2) {
			return(1);
		}
		else {
			return(0);
		}
	
	}
	
	
	
	public float get_transcript_to_sum_frag_counts (List<Integer> path) {
		if (! transcript_to_sum_frag_counts.containsKey(path)) {
			throw new RuntimeException("Error, path is not stored: " + path);
		}
		return(transcript_to_sum_frag_counts.get(path));
	}
	
	
	public static void main (String args[]) {
		
		// unit test
		BFLY_GLOBALS.VERBOSE_LEVEL = 20;
		
		
		// create transcripts:
		List<List<Integer>> transcripts = new ArrayList<List<Integer>>();
		
		//List<Integer>transcript_1 = new ArrayList<Integer>();
		Integer t1 [] = {1, 2, 3, 4, 5};
		List<Integer> transcript_1 = new ArrayList<Integer>(Arrays.asList(t1));
		
		Integer t2[] = {6, 2, 3, 4, 7};
		List<Integer> transcript_2 = new ArrayList<Integer>(Arrays.asList(t2));
		
		transcripts.add(transcript_1);
		transcripts.add(transcript_2);
		
		// set their lengths
		HashMap<List<Integer>, Integer> seqLengths = new HashMap<List<Integer>,Integer>();
		seqLengths.put(transcript_1, 1000);
		seqLengths.put(transcript_2, 1000);
		
		// add some reads in the form of pairpaths
		
		HashMap<List<Integer>, HashMap<PairPath,Integer>> transcript_to_pairpaths = new HashMap<List<Integer>, HashMap<PairPath,Integer>>();
		HashMap<PairPath,Integer> t1_to_read_support = new HashMap<PairPath,Integer>();
		HashMap<PairPath,Integer> t2_to_read_support = new HashMap<PairPath,Integer>();
		
		transcript_to_pairpaths.put(transcript_1, t1_to_read_support);
		transcript_to_pairpaths.put(transcript_2, t2_to_read_support);
		
		// p1 is shared between t1 and t2
		PairPath p1 = new PairPath();
		Integer p1_nodes [] = {2,3,4};
		p1.setPath1(new ArrayList<Integer>(Arrays.asList(p1_nodes)));
		t1_to_read_support.put(p1, 10);
		t2_to_read_support.put(p1, 10);
		
		
		
		// p2 is t1 only
		PairPath p2 = new PairPath();
		Integer p2_frag1_nodes [] = {1,2};
		Integer p2_frag2_nodes [] = {4,5};
		p2.setPath1(new ArrayList<Integer>(Arrays.asList(p2_frag1_nodes)));
		p2.setPath2(new ArrayList<Integer>(Arrays.asList(p2_frag2_nodes)));
		t1_to_read_support.put(p2, 8);
		
		// p3 is t2 only
		PairPath p3 = new PairPath();
		Integer p3_frag1_nodes [] = {6,2};
		Integer p3_frag2_nodes [] = {4,7};
		p3.setPath1(new ArrayList<Integer>(Arrays.asList(p3_frag1_nodes)));
		p3.setPath2(new ArrayList<Integer>(Arrays.asList(p3_frag2_nodes)));
		t2_to_read_support.put(p3, 2);
		
		PathExpressionComparator pc = new PathExpressionComparator(transcripts, transcript_to_pairpaths, seqLengths); 
		
		
		System.exit(0);
		
	}
	
	
		
}
	
	
	
