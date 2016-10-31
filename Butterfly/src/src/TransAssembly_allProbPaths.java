
//import java.util.Set;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistanceWoVer;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;
import gnu.getopt.*;

import java.io.BufferedReader;
//import java.io.FileNotFoundException;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.io.File;

import jaligner.Alignment;
//import org.jgrapht.alg.DijkstraShortestPath;

public class TransAssembly_allProbPaths {

	private static final boolean DEBUG = true;
	private static final SeqVertex ROOT = new SeqVertex(-1, "S",Integer.MAX_VALUE);
	private static final SeqVertex T_VERTEX = new SeqVertex(-2, "E");
	private static final int LINE_LEN = 60;
	private static int LAST_ID = -1;
	private static int LAST_REAL_ID = -1;
	private static int MAX_DEPTH = 0;

	private static double EDGE_THR = 0.05; // compares between each edge and its sister edges (u->v; vs all output of u, or all input of v)
	private static double FLOW_THR = 0.02;// compares between each edge and its flow of its vertices (u->v; vs all input of u, or all output of v)
	private static boolean NO_GRAPH_PRUNING = false; 
	
	private static final int COMP_AVG_COV_THR = 1;
	private static final int INITIAL_EDGE_ABS_THR = 0;


	private static int MIN_READ_SUPPORT_THR = 1;
	private static int MIN_OUTPUT_SEQ;

	// Paths Too Similar Settings
	private static int MAX_DIFFS_SAME_PATH = 2;
	private static int DIFFS_WINDOW_SIZE = 100;
	private static int MAX_FINAL_DIFFS_IN_WINDOW = 5;
	private static float MIN_PERCENT_IDENTITY_SAME_PATH = 98.0f;
	private static int MAX_INTERNAL_GAP_SAME_PATH = 10; 
	private static boolean NO_PATH_MERGING = false; // disables use of the above
	private static boolean NO_REMOVE_LOWER_RANKED_PATHS = false;
	private static boolean NO_GRAPH_COMPACTION = false; // also assumes no path merging

	
	// pasafly params
	private static int MAX_VERTEX_DISCOVERY_TIME_DIFF_ALLOW_COMPARE = 0; // used in pasafly-mode  //FIXME: should be zero, must investigate
	private static boolean FAST_PASA = false;
	private static int MAX_NUM_PATHS_PER_PASA_NODE = 10;
	
	// Path alignment settings
	//	private static int MAX_INTERNAL_GAP_LENGTH = 20; // minimum cassette exon size that might be skipped in an alt-splice variant.
	private static boolean SMITH_WATERMAN_ALIGN_FLAG = false;
	private static int MAX_SEQ_LEN_DP_ALIGN = 10000; 

	private static boolean MISO_OUTPUT = true;
	//private static boolean USE_PATH_ALIGNMENT = true;

	
	private static int MAX_PAIR_DISTANCE = 0;
	private static int PATH_REINFORCEMENT_DISTANCE_PERCENT = 25; 
	private static int PATH_REINFORCEMENT_DISTANCE = 0;
	private static int MAX_NUM_PATHS_PER_NODE_INIT = 100;
	private static int MAX_NUM_PATHS_PER_NODE_EXTEND = 25;
	
	// read sequence to graph mapping criteria
	private static int MAX_MM_ALLOWED = 0; // dynamically updated global (bad)
	private static int MAX_MM_ALLOWED_CAP = 0; // dynamically updated global (bad)
	private static double MAX_READ_SEQ_DIVERGENCE = 0.05;
	private static final double MAX_READ_LOCAL_SEQ_DIVERGENCE = 0.1; // avoid very bad locally aligned regions along the way.
	
	
	private static final int EXTREME_EDGE_FLOW_FACTOR = 200;

	// path extension alternative options
	private static final boolean USE_TRIPLETS = false; // do not use.  
	private static boolean ALL_POSSIBLE_PATHS = false; // most lenient form of path validation: all edge combinations allowed.
	private static boolean LENIENT_PATH_CHECKING = false; //  lenient: give benefit of doubt for connections that do not conflict

	// path reinforcement check options
	private static boolean ORIGINAL_PATH_EXTENSIONS = false; // examines paths from nodes to sinks
	private static int KMER_SIZE = 0;
	private static boolean GENERATE_FULL_SEQ_GRAPH = false;
	private static boolean GENERATE_MIDDLE_DOT_FILES = false;
	private static boolean COLLAPSE_SNPs = true;
	private static boolean TRIPLET_LOCKING = true;  //deprecated
	private static boolean EXTENDED_TRIPLET = true; // deprecated
	private static boolean TRIPLET_STRICT = false;  // deprecated
	private static boolean FRACTURE_UNRESOLVED_XSTRUCTURE = false;  // potentially useful - deciding to just annotate these in the fasta headers.
	private static boolean INFER_UNRESOLVED_XSTRUCTURE_PATHS = false; // harder problem than anticipated.  Turn to true statistical inference.
	
	private static boolean FIND_ALSO_DIFF_PATHS = false;

	private static boolean USE_DEGENERATE_CODE = false;
	private static String[] LETTERS = new String[]{"A","C","G","T"};

	private static PrintStream ERR_STREAM;
	private static boolean USE_STDERR = false;
	private static boolean LOG_STDERR = false;

	private static boolean NO_CLEANUP = false; // if set to true, removes input files

	private static boolean RUN_ALL_VS_ALL_FILTER = false;
	private static int ALL_VS_ALL_MAX_DP_LEN = 1000;

	
	private static HashMap<PairPath,ArrayList<String>> LONG_READ_PATH_MAP;  // PairPath => ArrayList(long_reads_names)
	private static HashMap<String, PairPath> LONG_READ_NAME_TO_PPath; // string => PairPath

	private static boolean cufflinksOpt = false; // minimum path set (Maria's code)
	private static boolean CUFF_NO_EXTEND = false; // enable to report only disjoint chains w/o extensions
	
	private static boolean EXPAND_LOOPS = true;  //FIXME: parameterize this after testing.
	
	
	private static boolean pasaFlyOpt = false; // PASA algorithm using pairpaths
	
	private static boolean pasaFlyUniqueOpt = false; // variation on PASA
	
	private static boolean ILLUSTRATE_FINAL_ASSEMBLIES = true;
	
	private static boolean MAKE_PE_SE = false; 
	
	private static boolean NO_EM_REDUCE = false;
	
	private static Float MIN_ISOFORM_PCT_LEN_OVERLAP = 30f;
	
	private static Integer GENE_COUNTER = 0;
	
	private static Float MIN_RELATIVE_ISOFORM_EXPRESSION = 5f;
	
	private static Integer READ_END_PATH_TRIM_LENGTH = 0;
	
	private static String FILE = ""; //  updated to 'file' value below.
	
	private static Boolean TREAT_PAIRS_AS_SINGLE = false;
	
	private static Integer PATH_COUNTER = 0;
	
	private static Boolean USE_DP_READ_TO_VERTEX_ALIGN = true;
	
	
	public static Comparator<SeqVertex> SeqVertexIDorderComparator = new Comparator<SeqVertex>() {

		public int compare (SeqVertex v_a, SeqVertex v_b) {

			// want most highly supported pairpaths to sort descendingly

			Integer v_a_id = v_a.getID();
			Integer v_b_id = v_b.getID();
			
			if (v_a_id < v_b_id) {
				return(-1);
			}
			else if (v_a_id > v_b_id) {
				return(1);
			}
			else {
				return(0);
			}
		}

	};

	
	
	private static Map<String, String> DEGENERATE_CODE = new HashMap<String, String>() {
		private static final long serialVersionUID = 1L;

		{ 
			put("AG","R");
			put("CT","Y");
			put("CG","S");
			put("AT","W");
			put("GT","K");
			put("AC","M");
			put("CGT","B");
			put("AGT","D");
			put("ACT","H");
			put("ACG","V");
			put("ACGT","N");
		}
	};

	private static Map<String, String> DEGENERATE_CODE_REV = new HashMap<String, String>() {
		private static final long serialVersionUID = 1L;

		{ 
			put("R","AG");
			put("Y","CT");
			put("S","GC");
			put("W","AT");
			put("K","GT");
			put("M","AC");
			put("B","CGT");
			put("D","AGT");
			put("H","ACT");
			put("V","ACG");
			put("N","ACGT");
		}
	};

	//private static Map<String, AlignmentStats> NUM_MATCHES_HASH;
	private static Map<String, AlignmentStats> NUM_MISMATCHES_HASH;


	public static void main(String[] args) throws Exception 
	{
		long totalNumReads = 0;

		/*
		for(int i = 0; i < args.length; i ++)
			System.out.println(args[i]);
		*/

		String file = "";
		boolean printUsage = false;
		LongOpt[] longopts = new LongOpt[100]; // big enough we don't have to keep incrementing it as our option list grows.
		longopts[0] = new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h');
		longopts[1] = new LongOpt("use-degenerate-code", LongOpt.OPTIONAL_ARGUMENT, null, 1); 
		longopts[2] = new LongOpt("dont-collapse-snps", LongOpt.OPTIONAL_ARGUMENT, null, 'S'); 
		longopts[3] = new LongOpt("generate-full-sequence-graphs", LongOpt.OPTIONAL_ARGUMENT, null, 'G');
		longopts[4] = new LongOpt("stderr", LongOpt.OPTIONAL_ARGUMENT, null, 2);

		StringBuffer sb = new StringBuffer(0);
		longopts[5] = new LongOpt("edge-thr", LongOpt.OPTIONAL_ARGUMENT, sb, 'E'); 
		longopts[6] = new LongOpt("flow-thr", LongOpt.OPTIONAL_ARGUMENT, sb, 'W'); 
		longopts[7] = new LongOpt("min_per_id_same_path", LongOpt.OPTIONAL_ARGUMENT, null, 3);
		longopts[8] = new LongOpt("max_number_of_paths_per_node_init", LongOpt.OPTIONAL_ARGUMENT, null, 4);
		longopts[9] = new LongOpt("min_per_align_same_path", LongOpt.OPTIONAL_ARGUMENT, null, 5);
		longopts[10] = new LongOpt("SW", LongOpt.NO_ARGUMENT, null, 6); //SMITH_WATERMAN_ALIGNMENT_FLAG

		longopts[11] = new LongOpt("all_possible_paths", LongOpt.NO_ARGUMENT, null, 7); // hidden option, testing only
		longopts[12] = new LongOpt("lenient_path_extension", LongOpt.NO_ARGUMENT, null, 8); // hidden for now

		longopts[13] = new LongOpt("path_reinforcement_distance", LongOpt.OPTIONAL_ARGUMENT, null, 9);

		longopts[14] = new LongOpt("original_path_extension", LongOpt.OPTIONAL_ARGUMENT, null, 10);
		longopts[15] = new LongOpt("ZIPPER", LongOpt.NO_ARGUMENT, null, 11); // hidden for now
		longopts[16] = new LongOpt("NO_MISO_OUTPUT", LongOpt.NO_ARGUMENT, null, 12); // hidden for now
		longopts[17] = new LongOpt("max_diffs_same_path", LongOpt.OPTIONAL_ARGUMENT, null, 13);
		longopts[18] = new LongOpt("max_internal_gap_same_path", LongOpt.OPTIONAL_ARGUMENT, null, 14);
		longopts[19] = new LongOpt("generate_intermediate_dot_files", LongOpt.NO_ARGUMENT, null, 15);
		longopts[20] = new LongOpt("triplet-lock", LongOpt.NO_ARGUMENT, null, 16); // deprecated
		longopts[21] = new LongOpt("max_seq_len_DP_align", LongOpt.OPTIONAL_ARGUMENT, null, 17);
		longopts[22] = new LongOpt("no_cleanup", LongOpt.NO_ARGUMENT, null, 18);
		longopts[23] = new LongOpt("log_stderr", LongOpt.NO_ARGUMENT, null, 19);
		longopts[24] = new LongOpt("__REDUCE", LongOpt.NO_ARGUMENT, null, 20); // deprecated
		longopts[25] = new LongOpt("diffs_window_size", LongOpt.OPTIONAL_ARGUMENT, null, 21);
		longopts[26] = new LongOpt("max_final_diffs_in_window", LongOpt.OPTIONAL_ARGUMENT, null, 22);
		longopts[27] = new LongOpt("CuffFly", LongOpt.NO_ARGUMENT, null, 23);
		longopts[28] = new LongOpt("no_path_merging", LongOpt.NO_ARGUMENT, null, 24);
		longopts[29] = new LongOpt("no_pruning", LongOpt.NO_ARGUMENT, null, 25);
		longopts[30] = new LongOpt("no_compaction", LongOpt.NO_ARGUMENT, null, 26);
		longopts[31] = new LongOpt("triplet_strict", LongOpt.NO_ARGUMENT, null, 27);
		longopts[32] = new LongOpt("extended_triplet", LongOpt.NO_ARGUMENT, null, 28);
		longopts[33] = new LongOpt("MAKE_PE_SE", LongOpt.NO_ARGUMENT, null, 29);
		longopts[34] = new LongOpt("cuff_no_extend", LongOpt.NO_ARGUMENT, null, 30);
		longopts[35] = new LongOpt("PasaFly", LongOpt.NO_ARGUMENT, null, 31);
		longopts[36] = new LongOpt("FAST_PASA", LongOpt.NO_ARGUMENT, null, 32);
		longopts[37] = new LongOpt("max_num_paths_per_pasa_node", LongOpt.OPTIONAL_ARGUMENT, null, 33);
		longopts[38] = new LongOpt("max_number_of_paths_per_node_extend", LongOpt.OPTIONAL_ARGUMENT, null, 34);
		longopts[39] = new LongOpt("PasaFlyUnique", LongOpt.NO_ARGUMENT, null, 35);
		longopts[40] = new LongOpt("NO_EM_REDUCE", LongOpt.NO_ARGUMENT, null, 36);
		longopts[41] = new LongOpt("MIN_PCT_DOM_ISO_EXPR", LongOpt.OPTIONAL_ARGUMENT, null, 37);
		longopts[42] = new LongOpt("READ_END_PATH_TRIM_LENGTH", LongOpt.OPTIONAL_ARGUMENT, null, 38);
		longopts[43] = new LongOpt("TREAT_PAIRS_AS_SINGLE", LongOpt.OPTIONAL_ARGUMENT, null, 39);
		longopts[44] = new LongOpt("no_remove_lower_ranked_paths", LongOpt.OPTIONAL_ARGUMENT, null, 40);
		longopts[45] = new LongOpt("NO_DP_READ_TO_VERTEX_ALIGN", LongOpt.NO_ARGUMENT, null, 41);
		longopts[46] = new LongOpt("MAX_READ_SEQ_DIVERGENCE", LongOpt.OPTIONAL_ARGUMENT, null, 42);
		
		Getopt g = new Getopt("TransAssembly", args, "L:F:N:C:V:SGDhO:R:",longopts);
		int c;


		if (MAX_SEQ_LEN_DP_ALIGN < ALL_VS_ALL_MAX_DP_LEN) {
			ALL_VS_ALL_MAX_DP_LEN = MAX_SEQ_LEN_DP_ALIGN;
		}

		while ((c = g.getopt()) != -1)
		{
			switch(c)
			{
			case 1:
				USE_DEGENERATE_CODE = true;
				break;
			case 2:
				USE_STDERR = true;
				break;
			case 3:
				MIN_PERCENT_IDENTITY_SAME_PATH = Float.parseFloat(g.getOptarg());
				break;
			case 4:
				MAX_NUM_PATHS_PER_NODE_INIT = Integer.parseInt(g.getOptarg());
				break;
			case 5:
				// no op
				break;
			case 6:
				SMITH_WATERMAN_ALIGN_FLAG = true;
				break;
			case 7:
				ALL_POSSIBLE_PATHS = true;
				break;
			case 8:
				LENIENT_PATH_CHECKING = true;
				break;
			case 9:
				PATH_REINFORCEMENT_DISTANCE = Integer.parseInt(g.getOptarg());
				break;
			case 10:
				ORIGINAL_PATH_EXTENSIONS = true;
				break;
			case 11:
				// available
				break;
			case 12:
				MISO_OUTPUT = false; // hidden option, that will output in MISO format
				break;
			case 13:
				MAX_DIFFS_SAME_PATH = Integer.parseInt(g.getOptarg());
				break;
			case 14:
				MAX_INTERNAL_GAP_SAME_PATH = Integer.parseInt(g.getOptarg());
				break;
			case 15:
				GENERATE_MIDDLE_DOT_FILES = true;
				break;
			
			case 16:
				TRIPLET_LOCKING = true;
				break;
		
			case 17:
				MAX_SEQ_LEN_DP_ALIGN = Integer.parseInt(g.getOptarg());
				break;
			case 18:
				NO_CLEANUP = true;
				break;
			case 19:
				LOG_STDERR = true;
				break;
			case 20:
				//RUN_ALL_VS_ALL_FILTER = true;
				break;
			case 21:
				DIFFS_WINDOW_SIZE = Integer.parseInt(g.getOptarg());
				break;
			case 22:
				MAX_FINAL_DIFFS_IN_WINDOW = Integer.parseInt(g.getOptarg());
				break;
			case 23:
				//min_path_set
				cufflinksOpt = true;
				break;
			case 24:
				// no path merging flag
				NO_PATH_MERGING = true;
				break;
			case 25:
				// turn off graph pruning (useful for testing purposes)
				NO_GRAPH_PRUNING = true;
				break;
			case 26:
				NO_GRAPH_COMPACTION = true;
				NO_GRAPH_PRUNING = true;
				break;
			
			
			case 27:
				TRIPLET_STRICT = true;
				TRIPLET_LOCKING = true;
				break;
			
			
			
			case 28:
				EXTENDED_TRIPLET = true;
				TRIPLET_LOCKING = true;
				break;
			
			
			case 29:
				MAKE_PE_SE = true;
				break;
				
			case 30:
				CUFF_NO_EXTEND = true;
				break;
				
			case 31:
				pasaFlyOpt = true;
				break;
				
			case 32:
				FAST_PASA = true;
				break;
				
			case 33:
				MAX_NUM_PATHS_PER_PASA_NODE = Integer.parseInt(g.getOptarg());
				break;
				
			case 34:
				MAX_NUM_PATHS_PER_NODE_EXTEND = Integer.parseInt(g.getOptarg());
				break;
				
			case 35:
				pasaFlyUniqueOpt = true;
				break;
				
			case 36:
				NO_EM_REDUCE = true;
				break;
				
			case 37:
				MIN_RELATIVE_ISOFORM_EXPRESSION = Float.parseFloat(g.getOptarg());
				break;
				
			case 38:
				READ_END_PATH_TRIM_LENGTH = Integer.parseInt(g.getOptarg());
				break;
				
			case 39:
				TREAT_PAIRS_AS_SINGLE = true;
				break;
				
				
			case 40:
				NO_REMOVE_LOWER_RANKED_PATHS = true;
				break;
				
			case 41:
				USE_DP_READ_TO_VERTEX_ALIGN = false;
				break;
				
			case 42:
				MAX_READ_SEQ_DIVERGENCE = Float.parseFloat(g.getOptarg());
				break;
				
				
			case 'S':
				COLLAPSE_SNPs = false;
				break;
			case 'G':
				GENERATE_FULL_SEQ_GRAPH = true;
				break;

			case 'h':
				printUsage = true;
				break;
			case 'L':
				MIN_OUTPUT_SEQ = Integer.parseInt(g.getOptarg());
				break;

			case 'F':
				MAX_PAIR_DISTANCE = Integer.parseInt(g.getOptarg());
				break;

			case 'N':
				totalNumReads = Long.parseLong(g.getOptarg());
				break;
			case 'V':
				BFLY_GLOBALS.VERBOSE_LEVEL = Integer.parseInt(g.getOptarg());
				break;
			case 'C':
				file = g.getOptarg();
				FILE = file;
				break;
			case 'D':
				FIND_ALSO_DIFF_PATHS = true;
				break;
			case 'O':
				PATH_REINFORCEMENT_DISTANCE_PERCENT = Integer.parseInt(g.getOptarg());
				break;
			case 'R':
				MIN_READ_SUPPORT_THR = Integer.parseInt(g.getOptarg());
				break;


			case 0:
				switch(Integer.parseInt(sb.toString()))
				{
				case 'E':
					// compares between each edge and its sister edges (u->v; vs all output of u, or all input of v)
					EDGE_THR = Double.parseDouble(g.getOptarg());
					break;

				case 'W':
					// compares between each edge and its flow of its vertices (u->v; vs all input of u, or all output of v)
					FLOW_THR = Double.parseDouble(g.getOptarg());
					break;
				}
				break;
			case '?':
				printUsage = true;
				break; 
				//
			default:
				printUsage = true;
			}
		}

		if (LOG_STDERR)
			ERR_STREAM = new PrintStream(new FileOutputStream(file + ".err"));

		debugMes("Started",10);


		debugMes("using Path alignment for path comparisons", 5);
		debugMes("combine paths if (identity=(numberOfMatches/shorterLen) > " + MIN_PERCENT_IDENTITY_SAME_PATH+"%" +
				" or if we have <= " + MAX_DIFFS_SAME_PATH+ " mismatches) "
				+ "and if we have internal gap lengths <= " + MAX_INTERNAL_GAP_SAME_PATH
				, 5); 


		int path_checking_opt_count = 0;
		if (LENIENT_PATH_CHECKING) {
			debugMes("Path extension mode: lenient.", 5);
			path_checking_opt_count++;
		}
		if (ORIGINAL_PATH_EXTENSIONS) {
			debugMes("Path extension mode: original path extension.", 5);
			path_checking_opt_count++;
		}
		if (ALL_POSSIBLE_PATHS) {
			debugMes("Path extension mode: all possible paths.", 5);
			path_checking_opt_count++;
		}

		if (path_checking_opt_count > 1) {
			System.err.println("Error, cannot enable more than one path checking option.");
			printUsage = true;
		}

		if (cufflinksOpt || pasaFlyOpt || pasaFlyUniqueOpt) {
			debugMes("CuffFly or PasaFly selected - running in SE mode to avoid uncertain alignments that break DAG transitivity", 5);
			MAKE_PE_SE = true;
		}
		
		

		printUsage = printUsage 
				|| file.equals("") 
				|| totalNumReads==0 
				|| MAX_PAIR_DISTANCE == 0 
				|| MIN_READ_SUPPORT_THR < 1;

		if (printUsage)
		{
			System.err.println("");
			System.err.println("########################################################################################");
			System.err.println("#");
			System.err.println("# Required:");
			System.err.println("#  -N  <int>     total number of reads or fragment pairs");
			System.err.println("#  -L  <int>     min length for an assembled sequence to be reported");
			System.err.println("#  -F  <int>     maximum fragment length (extreme dist between paired ends)");
			System.err.println("#  -C  <string>  prefix for component/reads file");
			System.err.println("#  ");
			System.err.println("#  ");
			System.err.println("#  ");
			System.err.println("# Optional:");

			System.err.println("#  ");
			System.err.println("# Graph compaction:");
			System.err.println("#  --edge-thr=<double>                sets the threshold for keeping the edge (u->v), compared to all *output* of u, or all *input* of v");
			System.err.println("#                                        (default: 0.05).");
			System.err.println("#  --flow-thr=<double>                sets the threshold for keeping the edge (u->v), compared to all *input* of u, or all *output* of v");
			System.err.println("#                                        (default: 0.02).");
			System.err.println("#  --no_pruning                       disable pruning of graph based on above thresholds.");
			System.err.println("#  --no_compaction                    do not compact the graph");
			System.err.println("#");
				
			System.err.println("# SNP handling modes:");
			System.err.println("#  --use-degenerate-code              use degenerate DNA code ");
			System.err.println("#                                        (default: don't use degenerate DNA code).");
			System.err.println("#  --dont-collapse-snps               don't collapse SNPs into a single letter ");
			System.err.println("#                                        (default: collapse SNPs into a single letter).");
			System.err.println("#  ");
			
			System.err.println("# Read-specific parameters:");
			System.err.println("#   --max_number_of_paths_per_node_init <int>      maximum number of unique pairpaths that can begin at a given node. (default: " + MAX_NUM_PATHS_PER_NODE_INIT + ")");
			System.err.println("#  --MAKE_PE_SE                                    convert split pair paths into separate single paths.");
			System.err.println("#  --NO_DP_READ_TO_VERTEX_ALIGN                    do not use DP alignment to align reads to vertex, use faster gap-free alignment");
			
			System.err.println("# Butterfly Path extension reinforcement requirements");
			System.err.println("      The following options are ordered by decreasing stringency.");
			System.err.println("#  --original_path_extension          examines paths from nodes to sinks, can be very slow");
			System.err.println("#   /compatible_path_extension/       *DEFAULT MODE* read (pair) must be compatible and contain defined minimum extension support for path reinforcement.");
			System.err.println("#  --lenient_path_extension           only the terminal node pair(v-u) require read support");
			System.err.println("#  --all_possible_paths               all edges are traversed, regardless of long-range read path support");
			System.err.println("#  --CuffFly                          cufflinks style assembly of minimum paths.");
			System.err.println("#  --cuff_no_extend                   do not extend max-matching pairpaths");
			
			System.err.println("#  -R <int>                           minimum read support threshold. Default: 2");
			System.err.println("#  -O <int>                           path reinforcement 'backwards overlap' percent of -F.  Default: (" + PATH_REINFORCEMENT_DISTANCE_PERCENT + ") Not used in --lenient_path_extension mode.");
			System.err.println("#      or");
			System.err.println("#  --path_reinforcement_distance=<int>  path reinforcement distance specified directly instead of computing it based on the -O value above.");
			//System.err.println("#  --triplet-lock                     lock in local (triplet) paths at nodes where read paths are supported. (increases stringency)");
			//System.err.println("#  --extended_triplet  (implies --triplet-lock) where read paths extend beyond a triplet, require path compatibility with extended read path prefix");
			//System.err.println("#  --triplet_strict    (implies --triplet-lock) break reconstructed paths at unsupported triplets");
			System.err.println("#  --max_number_of_paths_per_node_extend <int>    maximum number of paths that can be extended from a given node. (default: " + MAX_NUM_PATHS_PER_NODE_EXTEND + ")");
			System.err.println("#  --READ_END_PATH_TRIM_LENGTH <int>              min length of read terminus to extend into a graph node for it to be added to the pair path node sequence. (default: " + READ_END_PATH_TRIM_LENGTH + ")");
			System.err.println("#  --TREAT_PAIRS_AS_SINGLE              ignores pairing info.");
			System.err.println("#  ");
			System.err.println("#  ");
			
			
			System.err.println("# PasaFly related parameters:");
			System.err.println("#  --PasaFly                              use PASA-style pair path assembly");
			System.err.println("#  --FAST_PASA                            faster PASA by using a banded-DP strategy (experimental)");
			
			System.err.println("#  ");
			
			
			
			
			
			
			
			System.err.println("# Similar path reduction criteria:");
			System.err.println("#  --SW                                use Smith-Waterman local alignment mode (by default, uses Needleman-Wunsch global alignment)");
			System.err.println("#  --max_diffs_same_path=<int>         max allowed differences encountered between path sequences to combine them. (default: " + MAX_DIFFS_SAME_PATH + ")");
			System.err.println("#  --min_per_id_same_path=<float>      min percent identity for two paths to be merged into single paths (default: " + MIN_PERCENT_IDENTITY_SAME_PATH + ")");
			System.err.println("#  --max_internal_gap_same_path=<int>  maximum number of internal consecutive gap characters allowed for paths to be merged into single paths. (default: " + MAX_INTERNAL_GAP_SAME_PATH);
			System.err.println("#  --max_seq_len_DP_align=<int>        path regions to be aligned that are longer than this length use simple zipper alignment. (default: " + MAX_SEQ_LEN_DP_ALIGN + ")");
			System.err.println("#  --no_path_merging                   disable path comparisons and merging based on above settings.");
			//System.err.println("#  --no_remove_lower_ranked_paths      iteratively rank paths by unique read content, discard paths without unique support");
			
			System.err.println("#  --NO_EM_REDUCE                         run expectation maximization scheme to rank transcripts, remove lower ranking transcripts that yield no unique read content.");
			System.err.println("#  --MIN_PCT_DOM_ISO_EXPR=<float>      min percent expression of a dominantly expressed isoform for a gene. Default: " + MIN_RELATIVE_ISOFORM_EXPRESSION);
			System.err.println("#  ");

			System.err.println("# Misc: ");
			System.err.println("#  --max_number_of_paths_per_node     total number of paths allowed per node (default=" + MAX_NUM_PATHS_PER_NODE_INIT + ")");
			System.err.println("#  --generate-full-sequence-graphs    generate full sequence dot files");
			System.err.println("#                                        (default: generate dot files with start and end of each seq).");
			System.err.println("#  --generate_intermediate_dot_files  generate dot files after each step of butterfly");
			System.err.println("#                                        (default: generate only a the start and end dot files).");			
			System.err.println("#  --stderr                           prints the output to STDERR ");
			System.err.println("#  --log_stderr                       writes status info to fiule COMPONENT_PREFIX.err ");
			System.err.println("#  -V <int>                           verbosity level ");
			System.err.println("#                                        (default: 10 - progress of method + some stats)");
			System.err.println("#                                        (15 - like (10) + final paths to be added + additional loop info and dot files)");
			System.err.println("#                                        (20 - maximum verbosity)");
			System.err.println("#");
			System.err.println("########################################################################################");
			System.err.println("");
			System.exit(1);

		}
		
		if (USE_STDERR && BFLY_GLOBALS.VERBOSE_LEVEL < 5) {
			BFLY_GLOBALS.VERBOSE_LEVEL = 5; // default verbose level for stderr
		}
		
		
		// set calculated vars:
		if (PATH_REINFORCEMENT_DISTANCE > 0) {
			debugMes("path reinforcement distance set manually to: " + PATH_REINFORCEMENT_DISTANCE, 5);
		}
		else {

			PATH_REINFORCEMENT_DISTANCE = (int) (PATH_REINFORCEMENT_DISTANCE_PERCENT/100.0 * MAX_PAIR_DISTANCE);
			debugMes("path reinforcement distance computed based on " 
					+ PATH_REINFORCEMENT_DISTANCE_PERCENT + "% of max pair distance: "+ MAX_PAIR_DISTANCE
					+ " = " + PATH_REINFORCEMENT_DISTANCE + " bases", 5);
		}
		/* original version:
		if (PATH_REINFORCEMENT_DISTANCE == 0 && MAX_PAIR_DISTANCE > 50) {
			PATH_REINFORCEMENT_DISTANCE = MAX_PAIR_DISTANCE - 50; // Moran's original settings.
		}
		 */


		if (!COLLAPSE_SNPs && USE_DEGENERATE_CODE)
			USE_DEGENERATE_CODE = false;

		Vector<Integer> rootIDs = new Vector<Integer>();

		HashMap<Integer,Integer> outFlow = new HashMap<Integer, Integer>();
		HashMap<Integer,Integer> inFlow = new HashMap<Integer, Integer>();
		HashMap<Integer,String> kmers = new HashMap<Integer, String>();

		LONG_READ_PATH_MAP = new HashMap<PairPath, ArrayList<String>>();
		LONG_READ_NAME_TO_PPath = new HashMap<String, PairPath>();

		

		PrintStream pout_diff = null;
		PrintStream pout_all = new PrintStream(new FileOutputStream(file+".allProbPaths.fasta"));

		debugMes("SECTION\n================\nParsing de Bruijn graph\n======================\n", 5);
		debugMes("preProcessGraphFile: " + file + ".out", 10);
		preProcessGraphFile(file+".out",outFlow, inFlow, kmers);

		
		debugMes("SECTION\n==================\nbuildNewGraph\n========================\n", 5);
		debugMes("buildNewGraphFirstLetter: " + file + ".out", 10);
		DirectedSparseGraph<SeqVertex, SimpleEdge> graph = buildNewGraphUseKmers(file+".out",rootIDs,outFlow,inFlow,kmers); 

		SeqVertex.set_kmer_length(KMER_SIZE);
		SeqVertex.set_graph(graph);
		
		LAST_REAL_ID = LAST_ID;
		debugMes("Graph is built",10);

		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
			describeNodes(graph);
		}

		HashMap<String,Integer> originalGraphKmerToNodeID = new HashMap<String,Integer>();
		for (SeqVertex sv : graph.getVertices()) {
			String kmer = sv.getName();
			int id = sv.getID();
			originalGraphKmerToNodeID.put(kmer,  id);
			debugMes("ORIGINAL GRAPH NODE: " + kmer + " with ID: " + id, 20);
		}
		
		
		String[] tmpFile = file.split("/");
		String graphName = tmpFile[tmpFile.length-1];


		//boolean createMiddleDotFiles = false;
		boolean createMiddleDotFiles = GENERATE_MIDDLE_DOT_FILES;

		if (createMiddleDotFiles) 
			writeDotFile(graph,file + "_deBruijn.A.dot", graphName, false);
	
		
		// remember the original edge weights so we can relabel them later on in the final graph according to orig ids.
		HashMap<String,Double> original_edge_weights_using_orig_kmers = new HashMap<String,Double>();
		for (SimpleEdge e : graph.getEdges()) {
			String from_kmer = graph.getSource(e).getName();
			String to_kmer = graph.getDest(e).getName();
			
			original_edge_weights_using_orig_kmers.put(from_kmer + "_" + to_kmer, e.getWeight());
			
		}
		
		
		
		
		if (! NO_GRAPH_PRUNING) {
			debugMes("fixExtremeleyHighSingleEdges()", 1);
			fixExtremelyHighSingleEdges(graph,outFlow,inFlow);

			debugMes("removeLightEdges()", 1);
			removeLightEdges(graph);


			if (createMiddleDotFiles)
				writeDotFile(graph,file + "_removeLightEdges_init.B.dot",graphName, false);

		}
		
		
		if (! NO_GRAPH_COMPACTION) {
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
				debugMes("## Node descriptions before linear compaction:", 20);
				describeVertices(graph);
			}
			
			debugMes("compactLinearPaths()", 1);
			compactLinearPaths(graph);

			//removeShortOrphanNodes(graph, MIN_OUTPUT_SEQ);  // do this later
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
				debugMes("## Node descriptions after linear compaction:", 20);
				describeVertices(graph);
			}
			
			
			
			My_DFS dfs = new My_DFS(graph);
			dfs.runDFS2();
			
			
			if (createMiddleDotFiles)
				writeDotFile(graph,file + "_compactLinearPaths_init.C.dot",graphName, false);
		
			
			
		}
		
		
		/*
		
		if (! NO_GRAPH_PRUNING) { 
			boolean keep_pruning_graph = true;
			int prune_compact_round = 0;
			while (keep_pruning_graph) {

				prune_compact_round++;

				keep_pruning_graph = removeLightEdges(graph);

				compactLinearPaths(graph);

				if (createMiddleDotFiles)
					writeDotFile(graph,file + "_prune_compact_round_" + prune_compact_round + ".D.dot",graphName, false);

			}


		}
	   
		*/
		
		

		if (COLLAPSE_SNPs) {
			if (USE_DEGENERATE_CODE)
				removeSingleNtBubblesWithDegenerateCode(graph); // still needs updating //FIXME
			else
				removeSingleNtBubbles(graph);

			if (createMiddleDotFiles)
				writeDotFile(graph,file + "_SNPs_removed.H.dot",graphName, false);
			
			/*
			
			if (! NO_GRAPH_COMPACTION) {
				ompactLinearPaths(graph); 
			}

			if (createMiddleDotFiles)
				writeDotFile(graph,file + "_SNPs_removed_post-collapse.H.dot",graphName, false);
		
			*/
		
		}
		
		
		
		debugMes("SECTION\n====================\nRemoving small components.\n====================\n", 5);
		
		//remove small components
		calcSubComponentsStats(graph);

		if (graph.getVertexCount() == 0) {
			debugMes("Warning: graph pruned to nothingness", 1);
			//Runtime.getRuntime().exec("mv " + bfly_start_indicator_file + " " + bfly_end_indicator_file);
			System.exit(0);
		}

		if (createMiddleDotFiles)
			writeDotFile(graph,file + "_compactLinearPaths_removeSmallComp.D.dot",graphName, false);

		
		
		HashMap<Integer, LocInGraph> originalVerIDsMapping = getOriginalVerIDsMappingHash(graph);
		

		int numXstructs = countNumOfXstructures(graph);
		if (numXstructs>0)
			debugMes("number X structures = "+numXstructs,10);

		// Done Compacting graph.  

		DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis = new DijkstraDistance<SeqVertex, SimpleEdge>(graph, true);
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 19) {
			debugMes("\nSECTION\n=================\nNode descriptions before threading.\n===================\n", 5);
			describeNodes(graph);
		}
			
		
		
		
		debugMes("\nSECTION\n====================\nThreading reads through the graph\n=========================\n", 5);
		// maps individual reads to paths within the graph
		// readNameHash:  "actual read name" => Read object  (see Read class)
		HashMap<String, List<Read>> readNameHash = getReadStarts(graph,file+".reads",originalVerIDsMapping,rootIDs, originalGraphKmerToNodeID);

				
		

		debugMes("\nSECTION\n==================\nPairing up the reads into PairPaths\n===========================\n", 5);
		// note the logic for not doing the paired linking is handled under getReadStarts by just 
		// not using the /1 or /2 value.

		// pair up reads into PathPairs
		// combinedReadHash: start_vertex => (pair_path => count)
		HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash = getSuffStats_wPairs(graph,readNameHash,dijkstraDis);
		
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			debugMes("Printing Pair Paths  Before DAG Overlap Layout ------------------", 15);
			printPairPaths(combinedReadHash, "PairPaths@Init");
		}
		
		
		
		////////////////////////////////////////////////////////////////////////////////
		/// Moving from Collapsed de Bruijn Graph to an Overlap Graph 'seqvertex_graph' 
		////////////////////////////////////////////////////////////////////////////////
		
		
		debugMes("SECTION\n========  Create DAG from Overlap Layout ============\n\n", 5);

		DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph = new DirectedSparseGraph<SeqVertex, SimpleEdge>();

		HashMap<Integer, HashMap<PairPath, Integer>> seqvertex_combinedReadHash = create_DAG_from_OverlapLayout(seqvertex_graph, combinedReadHash, file, graphName, createMiddleDotFiles);

		

		My_DFS dfs = new My_DFS(seqvertex_graph);
		dfs.runDFS2();

		if (createMiddleDotFiles)
			writeDotFile(seqvertex_graph,file + "_vertex_DAG_postOverlapLayout.dot",graphName, false);


		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			debugMes("Printing Pair Paths ------------------", 15);
			printPairPaths(seqvertex_combinedReadHash, "PairPaths@PostOverlapLayout");
		}

		

		
		
		debugMes("SECTION\n======= Reorganize Read Pairings =========\n\n", 5);
		dijkstraDis = new DijkstraDistance<SeqVertex, SimpleEdge>(seqvertex_graph, true);
		seqvertex_combinedReadHash = reorganizeReadPairings(seqvertex_graph, seqvertex_combinedReadHash, dijkstraDis);
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15)
			printPairPaths(seqvertex_combinedReadHash, "PairPaths@AfterPairReorganization");
		
		
		//start working on one sub component at a time:
		// look for loops, try to solve them
		// if loops remain, move on to the next subComp.

		//removeShortOrphanNodes(graph, MIN_OUTPUT_SEQ);
		
		Set<Set<SeqVertex>> comps = divideIntoComponents(seqvertex_graph);   //**** IMPORTANT: THIS HAPPENS AFTER UNROLLING REPEATS AND BEFORE FINAL LOOP BREAKING

		debugMes("total number of components = "+comps.size(),10);
		int compID = -1;


		if (FIND_ALSO_DIFF_PATHS)
			pout_diff = new PrintStream(new FileOutputStream(file+"_diffPaths.fasta"));


		String[] pathName = file.split("/");


		if (seqvertex_graph.getVertexCount() == 0) {
			System.err.println("Graph is empty. Quitting.");
			System.exit(0);
		}


		int totalNumPaths = 0;
		int totalNumSuccComps = 0;
		
	
		// update the edge weights in the new graph
		for (SimpleEdge e : seqvertex_graph.getEdges()) {
			String from_kmer = seqvertex_graph.getSource(e).getLastKmer();
			String to_kmer = seqvertex_graph.getDest(e).getFirstKmer();
			
			debugMes("Searching for kmer set: " + from_kmer + " -> " + to_kmer, 15);
			
			Double orig_weight = original_edge_weights_using_orig_kmers.get(from_kmer + "_" + to_kmer);
			if (orig_weight != null)
				e.setWeight(orig_weight);
			
			
		}
		

		
		if (createMiddleDotFiles)
			writeDotFile(seqvertex_graph,file + "_final.Z.dot",graphName, false);
	
		
		////////////////////////////
		// Done with Graph Prep  ///
		////////////////////////////
		
		
		DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer = new DijkstraDistanceWoVer<SeqVertex, SimpleEdge>(seqvertex_graph);
 
		//initiation
		NUM_MISMATCHES_HASH = new HashMap<String, AlignmentStats>();  // cache alignment results for reuse.
        
        // --------------------------
		// Assemble transcript paths
        // --------------------------        
		
        for (Set<SeqVertex> comp : comps)
        {
        	compID++;
        	debugMes("\nSECTION\n============= Begin Assembly ===============\n\n", 5);
        	debugMes("Assembling subcomponent "+compID,10);
        	
        	totalNumSuccComps++;

        	// extract the reads that correspond to this subcomponent
        	HashMap<Integer,HashMap<PairPath,Integer>> componentReadHash = getComponentReads(compID, comp, seqvertex_combinedReadHash);

        	if (componentReadHash.isEmpty()) {
        		debugMes("No pairpaths stored for comp: " + comp + ", so skipping it.", 10);
        		continue;
        	}
        	
        	/*  old way... suffers from suboptimal setting of node depth
        	 * 
        	int count_pairpaths_removed = handleRemainingCyclicReads(componentReadHash, graph); 
        	
        	debugMes("Removed " + count_pairpaths_removed + " reads that appeared to retain complex cycles", 10);
        	*/
        	
        	if (BFLY_GLOBALS.VERBOSE_LEVEL >= 10) {
        		debugMes("####  Component Read Summary BEFORE PairPath-per-node Reduction", 10);
        		report_pairpath_counts(componentReadHash);
        	}
        	
        	reduce_to_max_paths_per_node(componentReadHash, TransAssembly_allProbPaths.MAX_NUM_PATHS_PER_NODE_INIT);
        	
        	if (BFLY_GLOBALS.VERBOSE_LEVEL >= 10) {
        		debugMes("####  Component Read Summary AFTER PairPath-per-node Reduction", 10);
        		report_pairpath_counts(componentReadHash);
        	}
        	
        	
        	// examine uncertainty of paths within the graph by looking at triplet support
        

        	debugMes("### Extracting triplets from reads.", 10);
        	HashMap<Integer, List<List<Integer>>> tripletMapper = extractTripletsFromReads(componentReadHash);
        	
        
        	HashMap<Integer,Boolean> xStructuresResolvedByTriplets = getXstructuresResolvedByTriplets(seqvertex_graph, comp, tripletMapper);

        	if (BFLY_GLOBALS.VERBOSE_LEVEL >= 10) {
        		// describe the locked down nodes
        		debugMes("\n### " + tripletMapper.size() + " nodes have locked-in triplet paths:", 10);
        		for (Integer central_node : tripletMapper.keySet()) {
        			debugMes("Triplet locks for: " + central_node + " : " + tripletMapper.get(central_node), 10);
        		}

        		

        	}

        	if (INFER_UNRESOLVED_XSTRUCTURE_PATHS) {
        		debugMes("## INFERRING UNRESOLVED X STRUCTURE PATHS ##", 10);
        		infer_best_triplets_across_unresolved_Xstructure(seqvertex_combinedReadHash, seqvertex_graph, xStructuresResolvedByTriplets, tripletMapper);
				
			}
        	
        	
        	
        	HashMap<List<Integer>, Pair<Integer>> FinalPaths_all = null;
        	
        	
        	
        	/*
        	if (false) { // just for debugging   //FIXME: make a debug parameter
        		
        		// just reconstruct paths based on the input reads
        		FinalPaths_all = reconstruct_paths_from_reads(graph, componentReadHash);

        	}
        	
        	else if (false) {
        		
        		// sort descendingly by length, and greedily assemble compatible paths.
        		
        		FinalPaths_all = reconstruct_paths_from_collapsed_reads(graph, componentReadHash);
        		
        		
        	}
        	else if (false) {
        		
        		examine_compatible_paths_debugging_only(graph, componentReadHash);
        		
        		System.exit(1);
        		
        		
        	}
        	
        	else if (false) {
        		
        		examine_out_of_order_depth_in_read_paths(graph, componentReadHash);
        		
        		System.exit(1);
        		
        	}
        	*/
        	
        	
        	
        		if (cufflinksOpt || pasaFlyOpt) {

        		
        		/*
        		// methods are very sensitive to out-of-order node depths in read paths:
        		int num_fractured_paths = handleRemainingCyclicReads(componentReadHash, graph);
        		
        		debugMes("Needed to fracture: " + num_fractured_paths + " pair paths due to out-of-order node depths", 10);
        		
        		*/
        		
        		if (cufflinksOpt) {

        			FinalPaths_all = cuffMinPaths(seqvertex_graph, componentReadHash,dijkstraDis);

        		}
        		else if (pasaFlyOpt) {

        			
                	debugMes("### Extracting complex path prefixes from reads.", 10);
                	HashMap<Integer, List<List<Integer>>> extendedTripletMapper = extractComplexPathPrefixesFromReads(componentReadHash);

        			
        			FinalPaths_all = pasafly(seqvertex_graph, componentReadHash,dijkstraDis, tripletMapper, extendedTripletMapper);


        		}
        	}
        	else if (pasaFlyUniqueOpt) {
        		
        		
            	debugMes("### Extracting complex path prefixes from reads.", 10);
            	HashMap<Integer, List<List<Integer>>> extendedTripletMapper = extractComplexPathPrefixesFromReads(componentReadHash);

        		
        		FinalPaths_all = pasaflyunique(seqvertex_graph, componentReadHash,dijkstraDis, tripletMapper, extendedTripletMapper);
        		
        	}
        	
        	else {
        		// Regular butterfly all probable paths mode:
        		
            	debugMes("### Extracting complex path prefixes from reads.", 10);
            	HashMap<Integer, List<List<Integer>>> extendedTripletMapper = extractComplexPathPrefixesFromReads(componentReadHash);

            	if (BFLY_GLOBALS.VERBOSE_LEVEL >= 16) {
            		debugMes("\n#### Extended triplets from reads: ", 16);
            		for (Integer term_node : extendedTripletMapper.keySet()) {
            			debugMes("Complex prefix paths for: " + term_node + " : " + extendedTripletMapper.get(term_node), 16);
            		}
            	}
            	
        		addSandT(seqvertex_graph,comp,componentReadHash);
        		
        		FinalPaths_all = butterfly(seqvertex_graph, comp, componentReadHash, totalNumReads, 
        				pout_all, dijkstraDis, dijkstraDisWoVer,
        				tripletMapper, extendedTripletMapper, xStructuresResolvedByTriplets);
        		//pathName = get_pathName_string(path, graph);
        	}
        	
        	
        	if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
        		for (List<Integer> path : FinalPaths_all.keySet()) {
        			debugMes("FinalPath@BeforeFiltering: " + path, 15);
        		}
        	}
        	
        	

        	// remove short paths
        	FinalPaths_all = remove_short_seqs(FinalPaths_all, seqvertex_graph);
        	
        	
        	if (FinalPaths_all.isEmpty()) {
        		debugMes("No paths to pursue. Continue...", 15);
        		continue;
        	}
        	
        	
        	int numXstructsResolved = countNumOfXstructuresResolved(seqvertex_graph,comp,FinalPaths_all);
        	if (numXstructs>0)
        		debugMes("number X structures resolved = "+numXstructsResolved + " / " + numXstructs,10);
        	
        	
        	// convert graph node IDs back to the original collapsed de Bruijn graph:
        	
        	HashMap<List<Integer>, Pair<Integer>> FinalPaths_all_orig_ids = new HashMap<List<Integer>, Pair<Integer>>();
        	for (List<Integer> final_path : FinalPaths_all.keySet()) {
        		List<Integer> revised_path_orig_ids = new ArrayList<Integer>();
        		for (Integer seq_vertex_id : final_path ) {
        			SeqVertex sv = SeqVertex.retrieveSeqVertexByID(seq_vertex_id);
        			
        			Integer orig_id = sv.getOrigButterflyID();
        			revised_path_orig_ids.add(orig_id);
        		}
        		FinalPaths_all_orig_ids.put(revised_path_orig_ids, FinalPaths_all.get(final_path));
        	}
        	

        	HashMap<List<Integer>,HashMap<PairPath,Integer>> finalPathsToContainedReads = assignCompatibleReadsToPaths(FinalPaths_all_orig_ids, combinedReadHash);

        	
        	if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
        		
        		// verbose dump of read support
        		
        		for (List<Integer> final_path : finalPathsToContainedReads.keySet()) {
        			HashMap<PairPath,Integer> contained_reads = finalPathsToContainedReads.get(final_path);
        			debugMes("PRELIM_FINAL_PATH:\n" + final_path + "\ncontains:", 20);
        			int sum_support = 0;
        			for (PairPath pp : contained_reads.keySet()) {
        				Integer read_support = contained_reads.get(pp);
        				debugMes(pp + "\tcount: " + read_support, 20);
        				sum_support += read_support;
        			}
        			debugMes("Total support: " + sum_support + "\n", 20);
        		}
        		
        	}
        	
        	// remove those paths that didn't have reads assigned:
        	{
        		Set<List<Integer>> paths_to_remove = new HashSet<List<Integer>>();
        		for (List<Integer> path : FinalPaths_all_orig_ids.keySet()) {
        			if (! finalPathsToContainedReads.containsKey(path)) {
        				debugMes("-removing final path that was not assigned read support: " + path, 10);
        				paths_to_remove.add(path);
        			}
        		}
        		for (List<Integer> path : paths_to_remove) {
        			FinalPaths_all_orig_ids.remove(path);
        		}
        	}
        	

        	
        	HashMap<List<Integer>, Pair<Integer>> filtered_paths_to_keep = new HashMap<List<Integer>,Pair<Integer>>();
        	
        	
        	/*
        	
        	if (! NO_REMOVE_LOWER_RANKED_PATHS) {


        		debugMes("SECTION\n======== Remove Lower Ranked Paths Without Unique Read Content ============\n\n", 5);

        		HashMap<List<Integer>, Pair<Integer>> lower_ranked_paths_removed = remove_lower_ranked_paths_without_unique_read_content(graph, FinalPaths_all, finalPathsToContainedReads);

        		filtered_paths_to_keep.putAll(lower_ranked_paths_removed);
        		
        	}

			*/
        	

        	
        	
        	//////////////////////////////////////
        	// Gene-level grouping of transcripts
        	//////////////////////////////////////
        	
        	
        	HashMap<List<Integer>,Integer> separate_gene_ids = group_paths_into_genes(FinalPaths_all_orig_ids, graph);
        
        	
        	//////////////////////////////////////
        	// Filtering out lower-quality paths
        	/////////////////////////////////////
        	
        	
        	
        	debugMes("Sep Gene IDs:" + separate_gene_ids, 10);
        	
        	if ( (! NO_EM_REDUCE) && FinalPaths_all.size() > 1) {

        		HashMap<List<Integer>, Pair<Integer>> EM_reduced_paths = run_EM_REDUCE(FinalPaths_all_orig_ids, graph, finalPathsToContainedReads, separate_gene_ids);

        		filtered_paths_to_keep.putAll(EM_reduced_paths);
        		
        	}
        	
        	// by default, running both lower ranking path removal and EM-reduction, and combining positively-filtered entries.
        	if (! filtered_paths_to_keep.isEmpty()) {
        		FinalPaths_all_orig_ids = filtered_paths_to_keep;
        	}
        	
        	if ( (! NO_PATH_MERGING)  && FinalPaths_all.size() > 1) {

        		// do CDHIT-like removal of highly similar but lesser supported paths.
        		debugMes("SECTION\n========= CD-HIT -like Removal of Too-Similar Sequences with Lesser Read Support =========\n\n", 5);

        		// alignment-based removal of lesser-supported paths that are too similar in sequence.
        		FinalPaths_all_orig_ids = reduce_cdhit_like(FinalPaths_all_orig_ids, graph, finalPathsToContainedReads);

        	}
        	
        	
        	
        	String component_name = pathName[pathName.length-1];


        	if (FinalPaths_all_orig_ids==null)
        		continue;


        	if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
        		for (List<Integer> path : FinalPaths_all_orig_ids.keySet()) {
        			debugMes("FinalPath@AfterFiltering: " + path, 15);
        		}
        	}
        	
        	
        	// get long read content:
        	HashMap<List<Integer>,ArrayList<String>> final_paths_to_long_read_content = new HashMap<List<Integer>,ArrayList<String>>();
        	if (! LONG_READ_PATH_MAP.isEmpty()) {
        		
        		assign_long_read_content_to_final_paths(FinalPaths_all_orig_ids, finalPathsToContainedReads,LONG_READ_PATH_MAP, final_paths_to_long_read_content);
        		
        	}
        	
        	
        	//-----------------------------
        	// Output the fasta sequences
        	//----------------------------
        	
        	printFinalPaths(FinalPaths_all_orig_ids, graph, pout_all, component_name, totalNumReads, 
        			final_paths_to_long_read_content, separate_gene_ids);
    		
        	
        	totalNumPaths += FinalPaths_all_orig_ids.size();

        	if ( BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
        		 debugMes("## ILLUSTRATING FINAL ASSEMBLIES", 20);
        		 illustrateFinalPaths(FinalPaths_all_orig_ids, finalPathsToContainedReads);
        	}
       
        	
        	
        }
       


		removeAllEdgesOfSandT(graph);

		pout_all.close();
		if (FIND_ALSO_DIFF_PATHS)
			pout_diff.close();


		debugMes("total number of paths reported = "+totalNumPaths+" from "+totalNumSuccComps +" components",1);

	

		debugMes("Done",10);
		
		if (LOG_STDERR)
			ERR_STREAM.close();
		
	}

	
	private static void assign_long_read_content_to_final_paths(
			HashMap<List<Integer>, Pair<Integer>> finalPaths_all,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> finalPathsToContainedReads,
			HashMap<PairPath, ArrayList<String>> LONG_READ_PATH_MAP_local,
			HashMap<List<Integer>, ArrayList<String>> final_paths_to_long_read_content) {


		for (List<Integer> path : finalPaths_all.keySet()) {

			if (finalPathsToContainedReads.get(path) == null) {
				//FIXME: why does this happen? very rare.
				continue;
			}
			
			for (PairPath pp : finalPathsToContainedReads.get(path).keySet()) {

				if (LONG_READ_PATH_MAP_local.containsKey(pp)) {

					if (! final_paths_to_long_read_content.containsKey(path)) {
						final_paths_to_long_read_content.put(path, new ArrayList<String>());

					}
					final_paths_to_long_read_content.get(path).addAll(LONG_READ_PATH_MAP_local.get(pp));

				}

			}



		}

		return;

	}






	private static HashMap<List<Integer>, Pair<Integer>> remove_short_seqs(
			HashMap<List<Integer>, Pair<Integer>> finalPaths_all,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		
		HashMap<List<Integer>, Pair<Integer>>  long_enough_paths = new HashMap<List<Integer>, Pair<Integer>>();
		
		for (List<Integer> path : finalPaths_all.keySet()) {
			
			String seq = getPathSeq(graph, path);
			
			if (seq.length() >=MIN_OUTPUT_SEQ) {
				
				// retain it:
				
				long_enough_paths.put(path, finalPaths_all.get(path));
				
			}
			
		}
		
		
		return(long_enough_paths);
		
	}


	public static HashMap<List<Integer>, Pair<Integer>> run_EM_REDUCE (HashMap<List<Integer>, Pair<Integer>> finalPaths_all, 
									  DirectedSparseGraph<SeqVertex, SimpleEdge> graph, HashMap<List<Integer>, 
									  HashMap<PairPath, Integer>> finalPathsToContainedReads, 
									  HashMap<List<Integer>,Integer> separate_gene_ids) {
		
		
		debugMes("SECTION\n====== ## BFLY_EM_REDUCE ## ==========\n\n", 5);
		
		
		List<List<Integer>> all_paths = new ArrayList<List<Integer>>(finalPaths_all.keySet());
		

		// need sequence lengths
		HashMap<List<Integer>,Integer> seqLengths = new HashMap<List<Integer>,Integer>();
		for (List<Integer> path : all_paths) {
			String seq = getPathSeq(graph,path);
			seqLengths.put(path,  seq.length());
		}


		//List<List<Integer>> retained_paths = new ArrayList<List<Integer>>();

		// sort paths by pair-path support descendingly
		PathExpressionComparator pc = new PathExpressionComparator(all_paths, finalPathsToContainedReads, seqLengths);



		Collections.sort(all_paths, pc);
		Collections.reverse(all_paths); // now descending according to read support.

		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			debugMes("Expression values for each candidate path:", 15);
			for (List<Integer> path : all_paths) {
				double expr = pc.get_expr(path);
				double sum_frag_counts = pc.get_transcript_to_sum_frag_counts(path);
				debugMes("Expr=" + expr + ", sum_exp_frags=" + sum_frag_counts + ", path: " + path, 15);

			}
		}


		// sort by expr, remove those w/ < 5% of expression of dominant isoform.
		List<List<Integer>> all_paths_min_rel_expr = remove_lesser_supported_paths_EM(all_paths, finalPathsToContainedReads, graph, pc, separate_gene_ids);
		
		
		// convert back to earlier-style data structure for compatibility
		HashMap<List<Integer>, Pair<Integer>> final_paths_map = new HashMap<List<Integer>, Pair<Integer>>();
		for (List<Integer> path : all_paths_min_rel_expr) {
			final_paths_map.put(path, finalPaths_all.get(path));
			debugMes("EM_REDUCE retaining: " + path, 15);
		}

		
		return (final_paths_map);



	}

	

	private static HashMap<Integer, HashMap<PairPath, Integer>> create_DAG_from_OverlapLayout(
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph, HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash, String dot_file_prefix, 
			String graphName, boolean createMiddleDotFiles) {
		
	
		
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();

		populate_pairpaths_and_readsupport(combinedReadHash, pairPaths, pairPathToReadSupport);

		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (PairPath pp : pairPaths) {

			paths.add(pp.getPath1());
			if (pp.hasSecondPath()) {
				paths.add(pp.getPath2());
			}
		}

		Collections.sort(paths, new Comparator<List<Integer>>() {

			public int compare(List<Integer> pathA, List<Integer> pathB) {

				if (pathA.size() < pathB.size()) {
					return(-1);
				}
				else if (pathA.size() > pathB.size()) {
					return(1);
				}
				else {
					return(0);
				}
			}



		});

		Collections.reverse(paths); // want descending by path l

		//////////////////////////////
		// remove the contained reads
		//////////////////////////////
		
		//contained_path_to_containers:  (key= the path contained, value = list of all other paths that fully contain it)
		HashMap<List<Integer>,List<List<Integer>>> contained_path_to_containers = new HashMap<List<Integer>,List<List<Integer>>>(); 
		List<List<Integer>> noncontained_paths = remove_containments(paths, contained_path_to_containers);
		
		debugMes("Noncontained paths: " + noncontained_paths, 15);
		
		
		//////////////////////////////
		// find dispersed repeats ////
		//////////////////////////////
		
		HashSet<Integer> dispersed_repeat_nodes = find_dispersed_repeat_nodes(noncontained_paths);
		
		
		////////////////////////////
		// build the overlap graph
		////////////////////////////
		
		
		// build a graph of compatible paths.
		List<Path> path_list = new ArrayList<Path>();
		for (List<Integer> p : noncontained_paths) {
			path_list.add(new Path(p));
		}
		
		HashMap<String,PathOverlap> pathMatches = new HashMap<String,PathOverlap>(); 
		
		DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph = construct_path_overlap_graph(path_list, pathMatches, 
				dispersed_repeat_nodes, dot_file_prefix, 
				 graphName,  createMiddleDotFiles);
		
		// draw the dot file for the path overlap graph:
		if (createMiddleDotFiles) 
			writeDotFile(path_overlap_graph, dot_file_prefix + "_POG.dot", graphName);
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			// output the path node listing
			for (Path p : path_overlap_graph.getVertices()) {
				debugMes("PathNodeDescription: " + p, 15);
			}
			
			
		}
		
		
		
		// add read pairing information to the graph:
		HashSet<SimplePathNodeEdge> pair_links = addPairPathsToOverlapGraph(path_overlap_graph, pairPathToReadSupport, contained_path_to_containers);
		
		// draw the dot file for the path overlap graph:
		if (createMiddleDotFiles) 
			writeDotFile(path_overlap_graph, dot_file_prefix + "_POG.PE_links_added.dot", graphName);

		
		//////////////////////////////
		// Breaking cycles
		/////////////////////////////
		
		int cycle_round = 0;
		
		boolean breaking_cycles = true;
		while (breaking_cycles) {
			
			cycle_round++;
			
			debugMes("// Breaking cycles in Path Overlap Graph (POG), Round: " + cycle_round, 10);
			
			breaking_cycles = break_cycles_in_path_overlap_graph(path_overlap_graph);
			
			if (createMiddleDotFiles) 
				writeDotFile(path_overlap_graph, dot_file_prefix + "_POG.cyclesRemoved.r" + cycle_round + ".dot", graphName);

		}
		
		// remove the pair_link edges before converting the overlap graph to a seq vertex graph, since 
		// edges in the overlap graph are intended to represent overlaps.
		for (SimplePathNodeEdge spne : pair_links) {
			if (path_overlap_graph.containsEdge(spne)) {
				path_overlap_graph.removeEdge(spne);
			}
		}
		
		//////////////////////////////////////////////
		// Convert the path DAG to a seq vertex DAG
		//////////////////////////////////////////////
		
		
		HashMap<Path,PathWithOrig> orig_path_to_updated_path = convert_path_DAG_to_SeqVertex_DAG(path_overlap_graph, 
				pathMatches, seqvertex_graph, dot_file_prefix, graphName,  createMiddleDotFiles);
		
		// note, path_overlap_graph includes non-contained paths
		//       pairPathToReadSupport contains all paths
		
		
		combinedReadHash = update_PairPaths_using_overlapDAG_refined_paths(orig_path_to_updated_path, pairPathToReadSupport, contained_path_to_containers);
		
		
		
		return(combinedReadHash);
		
		
	}


	private static HashSet<SimplePathNodeEdge> addPairPathsToOverlapGraph(
			DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph,
			Map<PairPath, Integer> pairPathToReadSupport,
			HashMap<List<Integer>, List<List<Integer>>> contained_path_to_containers) {
		
		// data structure conversion
		HashMap<List<Integer>,Path> list_to_path_hmap = new HashMap<List<Integer>,Path>();
		
		
		
		for (Path p : path_overlap_graph.getVertices()) {
			List<Integer> node_id_list = p.get_vertex_list();
			list_to_path_hmap.put(node_id_list, p);
		}
		
		
		HashSet<SimplePathNodeEdge> pair_link_edges = new HashSet<SimplePathNodeEdge> ();
		
		if (true)
			return(pair_link_edges);  // turning off PE //////////// DEBUG ////////////
		
		
		// add pair links
		for (PairPath pp : pairPathToReadSupport.keySet()) {
			if (pp.hasSecondPath()) {
				List<Integer> p1 = pp.getPath1();
				List<Integer> p2 = pp.getPath2();
				
				debugMes("# PE edges in overlap graph, targeting: " + p1 + " to " + p2, 15);
				
				
				List<Path> p1_path_list = new ArrayList<Path>();
				List<Path> p2_path_list = new ArrayList<Path>();
				
				if (list_to_path_hmap.containsKey(p1) && list_to_path_hmap.containsKey(p2)) {
					// add pairing edge between existing path nodes.
					p1_path_list.add(list_to_path_hmap.get(p1));
					p2_path_list.add(list_to_path_hmap.get(p2));
				}
				
				/*
				else if (list_to_path_hmap.containsKey(p1)) {
					
					p1_path_list.add(list_to_path_hmap.get(p1));
					
					// get containment list for p2, add edges from p1 -> p2
					
					List<List<Integer>> p2_contained_list = contained_path_to_containers.get(p2);
					for (List<Integer> p2_container : p2_contained_list) {
						Path p2_container_path = list_to_path_hmap.get(p2_container);
						p2_path_list.add(p2_container_path);
					}
					
				}
				else if (list_to_path_hmap.containsKey(p2)) {
					
					p2_path_list.add(list_to_path_hmap.get(p2));
					
					// get containment list for p1, and edges from p1-> p2
					
					List<List<Integer>> p1_contained_list = contained_path_to_containers.get(p1);
					for (List<Integer> p1_container : p1_contained_list) {
						Path p1_container_path = list_to_path_hmap.get(p1_container);
						p1_path_list.add(p1_container_path);
					}
					
				}
				*/

				int MAX_PAIR_LINKS = 3;
				if (p1_path_list.size() < MAX_PAIR_LINKS && p2_path_list.size() < MAX_PAIR_LINKS) {
					for (Path p1_path_node : p1_path_list) {
						
						for (Path p2_path_node : p2_path_list) {

							// add edge to graph if doesnt already exist:
							if (path_overlap_graph.findEdge(p1_path_node, p2_path_node) == null) {
								SimplePathNodeEdge spne = new SimplePathNodeEdge(1, p1_path_node.getPathNodeID(), p2_path_node.getPathNodeID());

								path_overlap_graph.addEdge(spne, p1_path_node, p2_path_node);	
								debugMes("-adding PE read edge between: " + p1_path_node.getPathNodeID() + " and " + p2_path_node.getPathNodeID(), 15);
							
								pair_link_edges.add(spne);
							}
						}
					}
				}
				
				
			} // endif pp.hasSecondPath()
		}
		
		
		return(pair_link_edges);
		
	}


	private static HashSet<Integer> find_dispersed_repeat_nodes(
			List<List<Integer>> paths) {
		
		HashMap<Integer,Integer> node_in_path_counter = new HashMap<Integer,Integer>();
		
		for (List<Integer> path : paths) {
			HashSet<Integer> node_found = new HashSet<Integer>();
			for (Integer node_id : path) {
				node_found.add(node_id);
			}
			
			Iterator<Integer> it = node_found.iterator();
			while (it.hasNext()) {
				Integer node_id = it.next();
				if (node_in_path_counter.containsKey(node_id)) {
					node_in_path_counter.put(node_id, node_in_path_counter.get(node_id)+1);
				}
				else {
					node_in_path_counter.put(node_id, 1);
				}
			}
			
		}
		
		
		
		
		
		List<Integer> node_ids = new ArrayList<Integer>(node_in_path_counter.keySet());
		
		final HashMap<Integer,Integer> node_counter = node_in_path_counter;
		
		Collections.sort(node_ids, new Comparator<Integer>() {
			public int compare (Integer a, Integer b) {
				if (node_counter.get(a) < node_counter.get(b)) {
					return(1);
				}
				else if (node_counter.get(a) > node_counter.get(b)) {
					return(-1);
					
				}
				else {
					return(0);
				}
			}
			
		});
		
		
		// pull out the repetitive ones
		
		int MIN_OCCURRENCE_REPEAT_NODE = 10;
		
		HashSet<Integer> repeat_nodes = new HashSet<Integer>();
		
		for (Integer node_id :  node_ids) {
			
			int repeat_count = node_counter.get(node_id);
			debugMes("Node[" + node_id + "] has repeat count: " + repeat_count, 15);
			
			if (repeat_count >= MIN_OCCURRENCE_REPEAT_NODE) {
				repeat_nodes.add(node_id);
			}
		}
		
		return(repeat_nodes);
		
	}


	private static HashMap<Integer, HashMap<PairPath, Integer>> update_PairPaths_using_overlapDAG_refined_paths(
			HashMap<Path, PathWithOrig> orig_path_to_updated_path,
			Map<PairPath, Integer> pairPathToReadSupport, 
			HashMap<List<Integer>, List<List<Integer>>> contained_path_to_containers) {
	
		
		// get the old-to-new listing in List<Integer> format for use with PairPath objects
		HashMap<List<Integer>,List<Integer>> old_to_new_path = new HashMap<List<Integer>,List<Integer>>();

		for (Path orig_path : orig_path_to_updated_path.keySet()) {
			List<Integer> orig_path_list = orig_path.get_vertex_list();

			List<Integer> updated_path_list = orig_path_to_updated_path.get(orig_path).getVertexList();

			old_to_new_path.put(orig_path_list, updated_path_list);

		}

		debugMes("Old-to-new-path mappings: " + old_to_new_path, 15);
		
		
		
		// get list of all old/new path pairs
		List<PathWithOrig> revised_paths = new ArrayList<PathWithOrig>(orig_path_to_updated_path.values());
		
		// now, create new pair paths based on updated mappings.
		HashMap<PairPath,Integer> updated_pairPaths = new HashMap<PairPath,Integer>();
		
		HashMap<PairPath,PairPath> old_pp_to_new_pp = new HashMap<PairPath,PairPath>();
		
		
		for (PairPath pp : pairPathToReadSupport.keySet()) {
			
			Integer read_support = pairPathToReadSupport.get(pp);
			
			debugMes("update_PairPaths_using_overlapDAG_refined_paths: orig_pp: " + pp + " has support: " + read_support, 20);
			
			
			PairPath new_pp;
			
			List<List<Integer>> p1_list = new ArrayList<List<Integer>>();
			
			List<Integer> p1 = pp.getPath1();
			if (old_to_new_path.containsKey(p1)) {
				p1_list.add(old_to_new_path.get(p1));
			}
			else {
				// might not be a unique path!! (eg. single original nodes now ending up in multiple places)
				p1_list = get_all_possible_updated_path_mappings(p1, revised_paths);
				
				debugMes("update_PairPaths_using_overlapDAG_refined_paths, p1: " + p1 + " mapped to: "  + p1_list, 20);
				
			}
			
			List<List<Integer>> p2_list = new ArrayList<List<Integer>>();
			
			if (pp.hasSecondPath()) {
				List<Integer> p2 = pp.getPath2();
				if (old_to_new_path.containsKey(p2)) {
					p2 = old_to_new_path.get(p2);
					p2_list.add(p2);
				}
				else {
					p2_list = get_all_possible_updated_path_mappings(p2, revised_paths);
				}
			
				// create new pair lists
				
				// restrict pair paths to those where each path maps uniquely
				if (p1_list.size() == 1 && p2_list.size() == 1) {
					List<Integer> p1_path = p1_list.get(0);
					List<Integer> p2_path = p2_list.get(0);
					new_pp = new PairPath(p1_path, p2_path);
					updated_pairPaths.put(new_pp, read_support);
					old_pp_to_new_pp.put(pp, new_pp);  // FIXME:  need to allow for multiple mappings here wrt long reads
				}
				else {
					// add each path separately if not already seen
					for (List<Integer> p1_path : p1_list) {
						if (! updated_pairPaths.containsKey(p1_path)) {
							new_pp = new PairPath(p1_path);
							updated_pairPaths.put(new_pp, 1);
						}
					}
					for (List<Integer> p2_path : p2_list) {
						if (! updated_pairPaths.containsKey(p2_path)) {
							new_pp = new PairPath(p2_path);
							updated_pairPaths.put(new_pp, 1);
						}
					}
					
					
				}
				/* orig
				for (List<Integer> p1_path : p1_list) {
					for (List<Integer> p2_path : p2_list) {
						new_pp = new PairPath(p1_path, p2_path);
						updated_pairPaths.put(new_pp, read_support);
						old_pp_to_new_pp.put(pp, new_pp);  // FIXME:  need to allow for multiple mappings here wrt long reads
						
					}
				}
				*/
			}
			else {
				// only individual paths
				for (List<Integer>p1_path : p1_list) {
					new_pp = new PairPath(p1_path);
					updated_pairPaths.put(new_pp, read_support);
					old_pp_to_new_pp.put(pp, new_pp); 
				}
			}
			
		}
		
		//update_long_read_path_mappings(old_pp_to_new_pp);
		
		
		HashMap<Integer, HashMap<PairPath, Integer>> new_combinedReadHash = construct_combinedReadhHash_from_PairPath_list(updated_pairPaths);
		
		return(new_combinedReadHash);
		
	}


	private static void update_long_read_path_mappings(
			HashMap<PairPath, PairPath> old_pp_to_new_pp) {
		
		
		debugMes("LONG_READ_PATH_MAP is:" + LONG_READ_PATH_MAP, 10);
		debugMes("LONG_READ_NAME_TO_PPath is : " + LONG_READ_NAME_TO_PPath, 10);
		
		
		
		HashMap<PairPath,ArrayList<String>> updated_LONG_READ_PATH_MAP = new HashMap<PairPath,ArrayList<String>>();  // PairPath => ArrayList(long_reads_names)
		HashMap<String, PairPath> updated_LONG_READ_NAME_TO_PPath = new HashMap<String,PairPath>(); // string => PairPath


		for (String long_read_name : LONG_READ_NAME_TO_PPath.keySet()) {
			PairPath pp = LONG_READ_NAME_TO_PPath.get(long_read_name);
			
			PairPath updated_pp = old_pp_to_new_pp.get(pp);
			
			updated_LONG_READ_NAME_TO_PPath.put(long_read_name, updated_pp);
			
			if (! updated_LONG_READ_PATH_MAP.containsKey(updated_pp)) {
				updated_LONG_READ_PATH_MAP.put(updated_pp, new ArrayList<String>());
			}
			updated_LONG_READ_PATH_MAP.get(updated_pp).add(long_read_name);
		
		}
		
		// replace old versions with updated versions.
		
		LONG_READ_PATH_MAP = updated_LONG_READ_PATH_MAP;
		LONG_READ_NAME_TO_PPath = updated_LONG_READ_NAME_TO_PPath;
		
		debugMes("LONG_READ_PATH_MAP updated to:" + updated_LONG_READ_PATH_MAP, 10);
		debugMes("LONG_READ_NAME_TO_PPath updated to : " + updated_LONG_READ_NAME_TO_PPath, 10);
		
		return;
	}


	private static List<Integer> update_path_mappings(List<Integer> p1,
			List<PathWithOrig> revised_paths) {
	
		PathWithOrig pwo_needs_updating = new PathWithOrig(p1);
		
		for (PathWithOrig pwo : revised_paths) {
			
			PathWithOrig updated_pwo = pwo_needs_updating.align_path_by_orig_id(pwo);
			if (updated_pwo != null) {
				return(updated_pwo.getVertexList());
			}
			
		}
		
		throw new RuntimeException("Unable to remap read: " + p1 + " given: " + revised_paths);
	}

	private static List<List<Integer>> get_all_possible_updated_path_mappings(
				List<Integer> p1,
				List<PathWithOrig> revised_paths) {
	
		List<List<Integer>> all_path_mappings = new ArrayList<List<Integer>>();
		
		PathWithOrig pwo_needs_updating = new PathWithOrig(p1);
		
		for (PathWithOrig pwo : revised_paths) {
			
			PathWithOrig updated_pwo = pwo_needs_updating.align_path_by_orig_id(pwo);
			if (updated_pwo != null) {
				List<Integer> updated_path = updated_pwo.getVertexList();
				if (! all_path_mappings.contains(updated_path)) {
					all_path_mappings.add(updated_path);
				}
			}
			
		}
		if (all_path_mappings.isEmpty()) {
		
			throw new RuntimeException("Unable to remap read: " + p1 + " given: " + revised_paths);
		}
		else {
			return(all_path_mappings);
		}
	}
	
	

	private static HashMap<Path,PathWithOrig> convert_path_DAG_to_SeqVertex_DAG(
			DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph,
			HashMap<String, PathOverlap> pathMatches,
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph, 
			String dot_file_prefix, 
			String graphName, 
			boolean createMiddleDotFiles) {
		
		
		debugMes("SECTION\n========  Convert Path-DAG to SeqVertex-DAG ============\n\n", 5);
		
		
	
		
		// init seqvertex graph to contain all nodes from expanded paths.
		
		HashMap<Path,List<SeqVertex>> orig_path_to_SeqVertex_list = new HashMap<Path,List<SeqVertex>>();
		HashMap<Path,PathWithOrig> orig_path_to_updated_path = new HashMap<Path,PathWithOrig>();
		
		for (Path p : path_overlap_graph.getVertices()) {
			List<Integer> node_id_list = p.get_vertex_list();
			List<SeqVertex> vertex_listing = new ArrayList<SeqVertex>();
			List<Integer> new_node_id_list = new ArrayList<Integer>();
			for (Integer node_id : node_id_list) {
				SeqVertex orig_vertex = SeqVertex.retrieveSeqVertexByID(node_id);
				Integer new_v_id = getNextID();
				SeqVertex new_v = new SeqVertex(new_v_id, orig_vertex);
				vertex_listing.add(new_v);
				seqvertex_graph.addVertex(new_v);
				
				new_node_id_list.add(new_v_id);
				
			}
			orig_path_to_SeqVertex_list.put(p, vertex_listing);
			
			PathWithOrig new_pwo = new PathWithOrig(p.getPathNodeID(), new_node_id_list, p.get_vertex_list());
			
			debugMes("prep_for_DAG_collapse: " + new_pwo, 15);
			
			orig_path_to_updated_path.put(p, new_pwo);
			
			// add edges between the vertices
			for (int i = 1; i < vertex_listing.size(); i++) {
				
				SeqVertex prev_v = vertex_listing.get(i-1);
				SeqVertex next_v = vertex_listing.get(i);
				
				SimpleEdge se = new SimpleEdge(1, prev_v.getID(), next_v.getID());
				seqvertex_graph.addEdge(se, prev_v, next_v);
				
			}
			
			
		}
		
		
		// do a DFS-based graph reconstruction starting from a root node.
		SeqVertex.set_graph(seqvertex_graph);
		
		
		HashSet<Path> visited = new HashSet<Path>();
		
		for (Path p : path_overlap_graph.getVertices()) {
			
			if (path_overlap_graph.getPredecessorCount(p) == 0) {
				// root node.
				DFS_add_path_to_graph(p, seqvertex_graph, path_overlap_graph, pathMatches, 
						orig_path_to_SeqVertex_list, visited);
			}
			
		}
		
		
		// before zippingUp
		if (createMiddleDotFiles)
			try {
				writeDotFile(seqvertex_graph, dot_file_prefix + "_before_zippingUpSeqVertexGraph.dot", graphName, false);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		

		if (graph_contains_loops(seqvertex_graph)) {
			throw new RuntimeException("Error, detected cycles in seqvertex_graph, so not a DAG as expected!");
		}

		List<SeqVertex> topo_sorted_vertices = TopologicalSort.topoSortSeqVerticesDAG(seqvertex_graph);
		

		
		
		// before zipping, after topo sort
		if (createMiddleDotFiles)
			try {
				writeDotFile(seqvertex_graph, dot_file_prefix + "_before_zippingUpSeqVertexGraph.TopoSort.dot", graphName, false);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		
		
		int zip_round = 0;
		
		int sum_merged = 1;
		
		while (sum_merged > 0) {
			sum_merged = 0;


			int count_zip_up_merged_in_round = 1;

			while (count_zip_up_merged_in_round > 0) {

				zip_round++;
				
				debugMes("\n\n## Round: " + zip_round + " Zipping up.", 10);
				
				if (graph_contains_loops(seqvertex_graph)) {
					throw new RuntimeException("Error, detected cycles in seqvertex_graph, so not a DAG as expected!");
				}

				init_replacement_vertices(seqvertex_graph);
				
				// ensure DAG 
				topo_sorted_vertices = TopologicalSort.topoSortSeqVerticesDAG(seqvertex_graph);


				count_zip_up_merged_in_round = zipper_collapse_DAG_zip_up(seqvertex_graph);
				sum_merged += count_zip_up_merged_in_round;

				debugMes("Zip up merged: " + count_zip_up_merged_in_round + " nodes.", 10);

				// draw the dot file for the path overlap graph:
				if (createMiddleDotFiles) {
					try {
						writeDotFile(seqvertex_graph, dot_file_prefix + "_zip_round_" + zip_round + "_zip_up.dot", graphName, false);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}


			int count_zip_down_merged_in_round = 1;

			while (count_zip_down_merged_in_round > 0) {

				zip_round++;
				debugMes("\n\n## Round: " + zip_round + " Zipping down.", 10);
				
				if (graph_contains_loops(seqvertex_graph)) {
					throw new RuntimeException("Error, detected cycles in seqvertex_graph, so not a DAG as expected!");
				}

				init_replacement_vertices(seqvertex_graph);
				
				// ensure DAG 
				topo_sorted_vertices = TopologicalSort.topoSortSeqVerticesDAG(seqvertex_graph);
				

				count_zip_down_merged_in_round = zipper_collapse_DAG_zip_down(seqvertex_graph);

				sum_merged += count_zip_down_merged_in_round;

				debugMes("Zip down merged: " + count_zip_down_merged_in_round + " nodes.", 10);

				// draw the dot file for the path overlap graph:
				if (createMiddleDotFiles) {
					try {
						writeDotFile(seqvertex_graph, dot_file_prefix + "_zip_round_" + zip_round + "_zip_down.dot", graphName, false);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}

			}

		}

			

		//  test again. :)
		if (graph_contains_loops(seqvertex_graph)) {
			throw new RuntimeException("Error, detected cycles in seqvertex_graph, so not a DAG as expected!");
		}

		// ensure DAG one last time
		topo_sorted_vertices = TopologicalSort.topoSortSeqVerticesDAG(seqvertex_graph);

		
		
		
		////////////////////////////////////////////////
		// update the paths based on their new vertices.
		////////////////////////////////////////////////
		
		// get old to new vertex id mapping
		HashMap<Integer,Integer> old_vertex_id_to_new_vertex_id = new HashMap<Integer,Integer>();
		for (SeqVertex v : topo_sorted_vertices) {
			
			Integer curr_vertex_id = v.getID();
			if (v.__tmp_compressed_vertices.size() > 0) {
				for (Integer old_vertex : v.__tmp_compressed_vertices) {
					old_vertex_id_to_new_vertex_id.put(old_vertex, curr_vertex_id);
					
					debugMes("Old_to_new_vertex_id_mapping: " + old_vertex + " => " + curr_vertex_id, 15);
				}
			}
			else {
				old_vertex_id_to_new_vertex_id.put(curr_vertex_id, curr_vertex_id);
				debugMes("Old_to_new_vertex_id_mapping: " + curr_vertex_id + " => " + curr_vertex_id + " (stays same)", 15);
			}
			
			
		}
		
		// update the old paths to the new paths
		for (PathWithOrig pwo : orig_path_to_updated_path.values()) {
			List<Integer> old_path = pwo.getVertexList();
			
			List<Integer> new_path = new ArrayList<Integer>();
			for (Integer id : old_path) {
				if (old_vertex_id_to_new_vertex_id.containsKey(id)) {
				Integer new_id = old_vertex_id_to_new_vertex_id.get(id);
				new_path.add(new_id);
				}
				else {
					throw new RuntimeException("Error, no new_id mapped from: " + id + ", in path: " + pwo);
				}
			}
			pwo.update_vertex_list(new_path);
		
		}
		
		
		return(orig_path_to_updated_path);
		
	
		
		
	}

	
	private static void init_replacement_vertices(
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph) {
		
		
		for (SeqVertex v : seqvertex_graph.getVertices()) {
			v.is_replacement_vertex = false;
		}
		
		return;
		
	}


	private static int zipper_collapse_DAG_zip_up(
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph) {

		
		int count_total_zip_merged = 0;
		
		/////////////////////////////
		// do bottom-up Zipping /////
		/////////////////////////////

		if (graph_contains_loops(seqvertex_graph)) {
			throw new RuntimeException("Error, detected cycles in seqvertex_graph, so not a DAG as expected!");	
		}

		List<SeqVertex> topo_sorted_vertices = TopologicalSort.topoSortSeqVerticesDAG(seqvertex_graph);

		Collections.reverse(topo_sorted_vertices);
		
		for (SeqVertex v : topo_sorted_vertices) {
			
			if (v.is_replacement_vertex) { continue; }
			
			if (! seqvertex_graph.containsVertex(v)) { continue; }
			
			
				count_total_zip_merged += zip_up(seqvertex_graph, v);
			
			
			
		}


		return(count_total_zip_merged);
		
		
	}

	
	
	private static int zipper_collapse_DAG_zip_down(
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph) {


		int count_total_zip_merged = 0;

		/////////////////////////////////////
		// do top-down zipping  /////////////
		/////////////////////////////////////


		if (graph_contains_loops(seqvertex_graph)) {
			throw new RuntimeException("Error, detected cycles in seqvertex_graph, so not a DAG as expected!");
		}

		List<SeqVertex> topo_sorted_vertices = TopologicalSort.topoSortSeqVerticesDAG(seqvertex_graph);




		for (SeqVertex v : topo_sorted_vertices) {

			if (v.is_replacement_vertex) { continue; }
			
			if (! seqvertex_graph.containsVertex(v)) { continue; }

			

			count_total_zip_merged += zip_down(seqvertex_graph, v);



		}


		return(count_total_zip_merged);


	}

	
	
	

	private static int zip_up(
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph,
			SeqVertex v) {
		
		
		
		List<SeqVertex> pred_list = new ArrayList<SeqVertex>(seqvertex_graph.getPredecessors(v));
		
		if (pred_list.size() <= 1) { return (0); } // must have multiple parents
		
		debugMes("## zip_up()", 15);
		
		// get list of parent nodes having the same original ID
		HashMap<Integer,HashSet<SeqVertex>> pred_orig_id_to_vertex_list = new HashMap<Integer,HashSet<SeqVertex>>();
		
		for (SeqVertex pred : pred_list) {
			
			if (pred.is_replacement_vertex) { return(0); } // delay to next round.
			
			if (! seqvertex_graph.containsVertex(pred)) { continue; }
			
			Integer orig_pred_id = pred.getOrigButterflyID();
			if (! pred_orig_id_to_vertex_list.containsKey(orig_pred_id)) {
				pred_orig_id_to_vertex_list.put(orig_pred_id, new HashSet<SeqVertex>());
			}
			pred_orig_id_to_vertex_list.get(orig_pred_id).add(pred);
			
		}

		
		int count_zip_merged = 0;
		
		
		for (HashSet<SeqVertex> pred_same_orig_id_set : pred_orig_id_to_vertex_list.values()) {
			
			if (pred_same_orig_id_set.size() == 1) { continue; } // need multiple parents for merging

			// merge them into a single node.
			
			

			count_zip_merged += attempt_zip_merge_SeqVertices(pred_same_orig_id_set, seqvertex_graph, "min");
		}
		
		return(count_zip_merged);
	}


	private static int zip_down (
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph,
			SeqVertex v) {
		
		
		
		List<SeqVertex> child_list = new ArrayList<SeqVertex>(seqvertex_graph.getSuccessors(v));
		
		if (child_list.size() <= 1) { return (0); } // must have multiple parents
		
		
		debugMes("##zip_down()", 15);
		
		
		// get list of children nodes having the same original ID
		HashMap<Integer,HashSet<SeqVertex>> child_orig_id_to_vertex_list = new HashMap<Integer,HashSet<SeqVertex>>();
		
		for (SeqVertex child : child_list) {
			
			if (child.is_replacement_vertex) { return(0); } // delay to next round
			
			if (! seqvertex_graph.containsVertex(child) ) { continue; }
			
			Integer orig_child_id = child.getOrigButterflyID();
			if (! child_orig_id_to_vertex_list.containsKey(orig_child_id)) {
				child_orig_id_to_vertex_list.put(orig_child_id, new HashSet<SeqVertex>());
			}
			child_orig_id_to_vertex_list.get(orig_child_id).add(child);
			
		}

		int count_zip_merged = 0;
		
		
		for (HashSet<SeqVertex> child_same_orig_id_set : child_orig_id_to_vertex_list.values()) {
			
			if (child_same_orig_id_set.size() == 1) { continue; } // need multiple parents for merging

			// merge them into a single node.
			
			

			count_zip_merged += attempt_zip_merge_SeqVertices(child_same_orig_id_set, seqvertex_graph, "max");
		}
		
		return(count_zip_merged);
	}
	
	
	
	private static int attempt_zip_merge_SeqVertices(HashSet<SeqVertex> pred_same_orig_id_set,
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph, String dir) {
	
		debugMes("attempt_zip_merge_SeqVertices(" + pred_same_orig_id_set + ")", 15);
		
		
		Integer replacement_vertex_id = getNextID();
		SeqVertex replacement_vertex_obj = null;
		
		// get list of all parents and all children of the merge-node-targets
		HashSet<SeqVertex> parent_vertices = new HashSet<SeqVertex>();
		HashSet<SeqVertex> child_vertices = new HashSet<SeqVertex>();

		HashSet<SimpleEdge> edges_to_delete = new HashSet<SimpleEdge>();
		
		// track depths, must ensure we keep the relative ordering in the DAG
		List<Integer> parent_depths = new ArrayList<Integer>();
		List<Integer> child_depths = new ArrayList<Integer>();
		List<Integer> target_depths = new ArrayList<Integer>();
		
		for (SeqVertex v : pred_same_orig_id_set) {
			
			int d = v.getNodeDepth();
			if (d < 0) {
				throw new RuntimeException("Error, seq vertex: " + v + " has negative depth setting");
			}
			
						
			target_depths.add(d);
			
			
			for (SeqVertex p : seqvertex_graph.getPredecessors(v)) {
				parent_vertices.add(p);
				
				parent_depths.add(p.getNodeDepth());
				
				// remove edge
				SimpleEdge se = seqvertex_graph.findEdge(p, v);
				edges_to_delete.add(se);
			}
			
			for (SeqVertex c: seqvertex_graph.getSuccessors(v)) {
				child_vertices.add(c);
				
				child_depths.add(c.getNodeDepth());
				
				// remove edge
				SimpleEdge se = seqvertex_graph.findEdge(v, c);
				edges_to_delete.add(se);
			}
			
			if (replacement_vertex_obj == null) {
				replacement_vertex_obj = new SeqVertex(replacement_vertex_id, v);
			}
			
			
			
		}
		
		
		if (parent_depths.size() > 0 && child_depths.size() > 0) {
		

			// ensure can merge and retain depth ordering:
			if ( ! (max_val(parent_depths) < min_val(child_depths) ) )
			{
				// cannot merge, since doing so would disrupt relative ordering of nodes
				return(0);
			}


		}
		
		
		
		// remove the graph edges:
		for (SimpleEdge se : edges_to_delete) {
			seqvertex_graph.removeEdge(se);
		}
		// remove the nodes themselves
		List<Integer> merged_vertex_ids = new ArrayList<Integer>();
		for (SeqVertex v : pred_same_orig_id_set) {
			merged_vertex_ids.add(v.getID());
			if (v.__tmp_compressed_vertices.size() > 0) {
				merged_vertex_ids.addAll(v.__tmp_compressed_vertices);
			}
			
			seqvertex_graph.removeVertex(v);
		}
		
		// add  new edges to parents
		for (SeqVertex p : parent_vertices) {
			SimpleEdge se = new SimpleEdge(1, p.getID(), replacement_vertex_obj.getID());
			seqvertex_graph.addEdge(se, p, replacement_vertex_obj);
		}
		
		// add new edges to children
		for (SeqVertex c : child_vertices) {
			SimpleEdge se = new SimpleEdge(1, replacement_vertex_obj.getID(), c.getID());
			seqvertex_graph.addEdge(se, replacement_vertex_obj, c);
		}
		
		Integer replacement_vertex_depth = (dir.equals("min")) ? min_val(target_depths) : max_val(target_depths);
		
		String zipDir = (dir.equals("min")) ? "Up" : "Down";
		
		replacement_vertex_obj.setDepth(replacement_vertex_depth);
		replacement_vertex_obj.setNodeDepth(replacement_vertex_depth);
		
		replacement_vertex_obj.is_replacement_vertex = true;
		
		// fix local environment for this round
		for (SeqVertex p : parent_vertices) {
			p.is_replacement_vertex = true;
		}
		for (SeqVertex c : child_vertices) {
			c.is_replacement_vertex = true;
		}
		
		replacement_vertex_obj.__tmp_compressed_vertices.addAll(merged_vertex_ids);
		
		debugMes(zipDir + "ZipMerging nodes: " + pred_same_orig_id_set + " to " + replacement_vertex_obj, 15);
		
		
		int count_merged = pred_same_orig_id_set.size();
		return(count_merged);
		
	}


	private static int max_val(List<Integer> vals) {
		
		
		Integer max_val = null;
		
		for (Integer val : vals) {
			if (max_val == null || val > max_val) {
				max_val = val;
			}
		}
		
		return(max_val);
		
	}

	
	private static int min_val(List<Integer> vals) {
		
		
		Integer min_val = null;
		
		for (Integer val : vals) {
			if (min_val == null || val < min_val) {
				min_val = val;
			}
		}
		
		return(min_val);
		
	}


	private static void DFS_add_path_to_graph(Path p,
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph,
			DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph, 
			HashMap<String, PathOverlap> pathMatches, 
			HashMap<Path, List<SeqVertex>> orig_path_to_SeqVertex_list, 
			HashSet<Path> visited) {
	
		
		if (visited.contains(p)) {
			// already done
			return;
		}
		
		debugMes("\nDFS_path_to_graph: targeting: " + p, 15);
		
		visited.add(p);
		
		//////////////////////////////////////////////////////////////////////////////////
		// Phase 1.  find candidate adjacent paths for use in labeling nodes in this path.
		
		
		List<Path> adjacent_untraversed_pathnodes = new ArrayList<Path>(); // for later, deciding next DFS entries
		
	
		
		for (Path succ : path_overlap_graph.getSuccessors(p)) {
			
			String pair_token = get_path_compare_token(p, succ);
			PathOverlap po = pathMatches.get(pair_token);
			
			int match_len = po.match_length;
			
			// draw edge between curr last node and next node in the successor path
			List<SeqVertex> curr_vertex_list = orig_path_to_SeqVertex_list.get(p);
			List<SeqVertex> succ_vertex_list = orig_path_to_SeqVertex_list.get(succ);
			
			
			boolean connect_all_matching_positions = true;
			if (connect_all_matching_positions) {


				for (int i = curr_vertex_list.size() - match_len, j = 0;
						i < curr_vertex_list.size() && j < match_len;
						i++,j++) {


					SeqVertex curr_vertex = curr_vertex_list.get(i);
					SeqVertex succ_vertex = succ_vertex_list.get(j+1);

					SimpleEdge se = new SimpleEdge(1, curr_vertex.getID(), succ_vertex.getID());
					seqvertex_graph.addEdge(se, curr_vertex, succ_vertex);


				}
			}
			else {
				
				// just the last one
				
				SeqVertex curr_vertex = curr_vertex_list.get(curr_vertex_list.size()-1);
				SeqVertex succ_vertex = succ_vertex_list.get(match_len); // linking up the prev to next+1

				SimpleEdge se = new SimpleEdge(1, curr_vertex.getID(), succ_vertex.getID());
				seqvertex_graph.addEdge(se, curr_vertex, succ_vertex);

				
				
			}

			
			DFS_add_path_to_graph(succ,
					seqvertex_graph,
					path_overlap_graph, 
					pathMatches, 
					orig_path_to_SeqVertex_list, 
					visited);
			
		}
		
		
		
		
		
		return;
		
	}
	
	

	/* orig code
	
	private static void DFS_add_path_to_graph(Path p,
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph,
			DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph,
			HashMap<Path, PathWithOrig> orig_path_to_updated_path, 
			HashMap<String, PathOverlap> pathMatches) {
	
		
		debugMes("\nDFS_path_to_graph: targeting: " + p, 15);
		
		//////////////////////////////////////////////////////////////////////////////////
		// Phase 1.  find candidate adjacent paths for use in labeling nodes in this path.
		
		
		List<Path> adjacent_untraversed_pathnodes = new ArrayList<Path>(); // for later, deciding next DFS entries
		
		// find a predecessor or successor that has the greatest overlap
		// and is already part of the new graph
		
		PathWithOrig best_predecessor_path = null;
		PathOverlap best_predecessor_overlap = null;
		
		
		
		for (Path pred : path_overlap_graph.getPredecessors(p)) {
			
			String pair_token = get_path_compare_token(pred, p);
			PathOverlap po = pathMatches.get(pair_token);
			
			if (orig_path_to_updated_path.containsKey(pred)) {
				// candidate for use as template for node assignment.
				
				if (best_predecessor_overlap == null || best_predecessor_overlap.match_score < po.match_score) {
					best_predecessor_overlap = po;
					best_predecessor_path = orig_path_to_updated_path.get(pred);
				}
			}
			else {
				pred._tmp_score = po.match_score;
				adjacent_untraversed_pathnodes.add(pred);
			}
		}
	
		
		PathWithOrig best_successor_path = null;
		PathOverlap best_successor_overlap = null;
		
		for (Path succ : path_overlap_graph.getSuccessors(p)) {
			
			String pair_token = get_path_compare_token(p, succ);
			PathOverlap po = pathMatches.get(pair_token);
			
			if (orig_path_to_updated_path.containsKey(succ)) {
				// candidate for use as template for node assignment
				
				if (best_successor_overlap == null || best_successor_overlap.match_score < po.match_score) {
					best_successor_overlap = po;
					best_successor_path = orig_path_to_updated_path.get(succ);
				}
			}
			else {
				succ._tmp_score = po.match_score;
				adjacent_untraversed_pathnodes.add(succ);
			}
		}
		
		///////////////////////////////////////////////////////////////////////
		// Phase 2: Refine labeling
		
		// init the new path:
		List<Integer> new_path = new ArrayList<Integer>();
		for (Integer i : p.get_vertex_list()) {
			new_path.add(-1);
		}
		
		PathWithOrig new_pwo = new PathWithOrig(p.getPathNodeID(), new_path, p.get_vertex_list());
		
		if (best_predecessor_path == null && best_successor_path == null) {
			
			debugMes("-dfs_msg: no best predecessor or successor path, so adding orig path from scratch.", 15);
			
			// start building the graph here.
			List<Integer> updated_path = new ArrayList<Integer>();
			SeqVertex prev_vertex = null;
			for (Integer orig_node_id : new_pwo.getOrigVertexList()) {
				SeqVertex v = SeqVertex.retrieveSeqVertexByID(orig_node_id);
				Integer next_v_id = getNextID();
				SeqVertex new_v = new SeqVertex(next_v_id,v);
				seqvertex_graph.addVertex(new_v);
				updated_path.add(next_v_id);
				if (prev_vertex != null) {
					// add new edge
					SimpleEdge se = new SimpleEdge(1, prev_vertex.getID(), next_v_id);
					seqvertex_graph.addEdge(se, prev_vertex, new_v);
				}
				prev_vertex = new_v;
				
			}
			new_pwo.update_vertex_list(updated_path);
			
			
		}
		else {
			
			
			List<Integer> updated_path = new_pwo.getVertexList(); // original vertex list
			
			
			
			// update nodes based on best matching predecessor
			
			if (best_predecessor_path != null) {
				debugMes("-dfs_msg: updating path " + p + " based on best predecessor: " + best_predecessor_path, 15);
				List<Integer> predecessor_node_ids = best_predecessor_path.getVertexList();
				for (int i = 0, j = predecessor_node_ids.size() - best_predecessor_overlap.match_length; 
						i < best_predecessor_overlap.match_length && j < predecessor_node_ids.size(); 
						i++, j++) {
					if (updated_path.get(i) != -1 && updated_path.get(i) != predecessor_node_ids.get(j)) {
						throw new RuntimeException("conflict in path assignments: " + updated_path + ", " + best_predecessor_path);
					}
					updated_path.set(i, predecessor_node_ids.get(j));
				}
			}
			
			
			// update nodes based on best matching successor

			
			if (best_successor_path != null) {
				debugMes("-dfs_msg: updating path " + p + " based on best successor: " + best_successor_path, 15);
				List<Integer> successor_node_ids = best_successor_path.getVertexList();
				for (int i = updated_path.size() - best_successor_overlap.match_length, j = 0;
						j < best_successor_overlap.match_length && i < updated_path.size();
						i++, j++) {

					if (updated_path.get(i) != -1 && updated_path.get(i) != successor_node_ids.get(j)) {
						throw new RuntimeException("conflict in path assignments: " + updated_path + "," + best_successor_path);
					}
					updated_path.set(i, successor_node_ids.get(j));
				}
			}
		   
			
			
			// add new nodes and edges for those that are path-specific here.
			for (int i = 0; i < updated_path.size(); i++) {
				List<Integer> orig_path = new_pwo.getOrigVertexList();
				
				if (updated_path.get(i) == -1) {
					// need new node:
					SeqVertex orig_v = SeqVertex.retrieveSeqVertexByID(orig_path.get(i));
					Integer new_node_id = getNextID();
					SeqVertex new_v = new SeqVertex(new_node_id, orig_v);
					seqvertex_graph.addVertex(new_v);
					updated_path.set(i, new_node_id);
				}
			}
			// ensure edges exist among this node set:
			for (int i = 1; i < updated_path.size(); i++) {
				SeqVertex prev_vert = SeqVertex.retrieveSeqVertexByID(updated_path.get(i-1));
				SeqVertex curr_vert = SeqVertex.retrieveSeqVertexByID(updated_path.get(i));
				SimpleEdge se = seqvertex_graph.findEdge(prev_vert, curr_vert);
				if (se == null) {
					se = new SimpleEdge(1, prev_vert.getID(), curr_vert.getID());
					seqvertex_graph.addEdge(se, prev_vert, curr_vert);
				}
				
			}
			
			
		}
		
		orig_path_to_updated_path.put(p, new_pwo);
		
		debugMes("-dfs_msg: newly added path is: " + new_pwo, 15);
		
		////////////////////////////////////////////////////////////////
		// phase 3:  DFS to next best overlapping adjacent edge.
		
		// get list of all edges not yet traversed
		// sort by match score
		// DFS them in order of match score
		
		Collections.sort(adjacent_untraversed_pathnodes, new Comparator<Path>() {
			
			public int compare (Path a, Path b) {
				if (a._tmp_score < b._tmp_score) {
					return(1);
				}
				else if (a._tmp_score > b._tmp_score) {
					return(-1);
				}
				else {
					return(0);
				}
			}
			
		});

		for (Path next_p : adjacent_untraversed_pathnodes) {
			if (! orig_path_to_updated_path.containsKey(next_p)) {
				DFS_add_path_to_graph(next_p, seqvertex_graph, path_overlap_graph, orig_path_to_updated_path, pathMatches);
			}
		}
		
		
		
		
		return;
		
	}

	*/


	private static boolean break_cycles_in_path_overlap_graph(
			DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph) {
	
		DijkstraShortestPath<Path, SimplePathNodeEdge> dp = new DijkstraShortestPath<Path, SimplePathNodeEdge>(path_overlap_graph);

		Set<Set<SimplePathNodeEdge>> curLoops = new HashSet<Set<SimplePathNodeEdge>>();

		// find all loops in the graph by seeing if, given edge v->v2, there is a path from v2 back to v

		for (Path p : path_overlap_graph.getVertices()) {
			for (Path s : path_overlap_graph.getSuccessors(p))
			{
				if (dp.getDistance(s, p)!=null) // there is a connection between p->s->....->p
				{
					//path has all edges from v to itself thru v2
					List<SimplePathNodeEdge> loopPath = dp.getPath(s, p);

					// v2 is successor of v, so let's just add the v->v2 edge too, complete the full loop.
					loopPath.add(0, path_overlap_graph.findEdge(p, s)); 

				
					// Collect the loop edge set.
					Set<SimplePathNodeEdge> loopPath_set = new HashSet<SimplePathNodeEdge>(loopPath);

					if (!curLoops.contains(loopPath_set))
					{
						curLoops.add(loopPath_set);
						debugMes("Found loop: " + loopPath_set, 15);
						
					}
				}
			}

		}
		
		if (curLoops.isEmpty())
			return false; // no cycles to break

		
		////////////////////////
		// process found loops
		///////////////////////


		
		
		Set<SimplePathNodeEdge> allRelevantEdges = new HashSet<SimplePathNodeEdge>();

		for (Set<SimplePathNodeEdge> loopPath_set : curLoops)
			for (SimplePathNodeEdge e : loopPath_set)
			{
				e.increaseNumOfLoopsBy1();
				allRelevantEdges.add(e);
			}


		////////////////////////////////////////
		//  break complex loops
		////////////////////////////////////////
		
		boolean res = false;

		if (!allRelevantEdges.isEmpty()){

			Comparator<SimplePathNodeEdge> numLoopsComparator = new NumPathNodeLoopsEdgeComparator();

			PriorityQueue<SimplePathNodeEdge> edgesQ = new PriorityQueue<SimplePathNodeEdge>(allRelevantEdges.size(), numLoopsComparator);

			edgesQ.addAll(allRelevantEdges);

			//while there are still loops
			// find the next edge that can be removed to reduce the number of loops
			// updated queue: remove all edges, and update their loop content

			SimplePathNodeEdge nextEtoRemove;
			while ( (!curLoops.isEmpty()) && (! edgesQ.isEmpty()) ) 
			{

				//FIXME: there was a situation where curLoops was not empty, 
				//       but edgesQ was, so I added edgesQ to the while condition. Investigate why this might happen.
				//       In this case, a node was involved in a self loop and a double-loop.

				nextEtoRemove = edgesQ.poll();

				if (path_overlap_graph.getSource(nextEtoRemove) == null 
						|| path_overlap_graph.getDest(nextEtoRemove) == null 
						|| nextEtoRemove.getNumOfLoopsInvolved() <= 0) { 
					continue; 
				}

				debugMes("removing the edge " + path_overlap_graph.getSource(nextEtoRemove).getPathNodeID() + "->" + 
						path_overlap_graph.getDest(nextEtoRemove).getPathNodeID() + " that appears in "
						+nextEtoRemove.getNumOfLoopsInvolved() + " loops",15);

				
				// remove the loops that have this edge from curLoops
				Set<Set<SimplePathNodeEdge>> removeLoops = new HashSet<Set<SimplePathNodeEdge>>();
				for (Set<SimplePathNodeEdge> loopPath_set : curLoops)
					if (loopPath_set.contains(nextEtoRemove))
					{
						debugMes("the loop "+ loopPath_set+" is now solved",15);
						removeLoops.add(loopPath_set);

						// update the number of loops involved in each edge
						for (SimplePathNodeEdge e : loopPath_set)
							e.decreaseNumOfLoopsBy1();
					}

				for (Set<SimplePathNodeEdge> loopPath_set : removeLoops)
					curLoops.remove(loopPath_set);


				//update the queue. remove all, and insert again if numLoops>0.
				SimplePathNodeEdge[] relEdges = (SimplePathNodeEdge[]) edgesQ.toArray(new SimplePathNodeEdge[0]);
				edgesQ.clear();
				for (SimplePathNodeEdge otherE : relEdges)
					if (otherE.getNumOfLoopsInvolved()>0)
						edgesQ.add(otherE);

				// remove this edge
				path_overlap_graph.removeEdge(nextEtoRemove);
				res = true;
			}
		}
		return res;

	
	
	}




	private static void writeDotFile(
			DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph,
			String output_filename, String graphName) {

		PrintStream p;
		try {
			p = new PrintStream(new FileOutputStream(output_filename));
			
			p.println("digraph G {");
			Path toVertex;

			//for each edge decide it's color
			for (Path vertex : path_overlap_graph.getVertices())
			{ //go over all vertices

				String verDesc = ""+vertex.getPathNodeID() +" [label=\"" + vertex.getPathNodeID() + "\"]";
				p.println("\t" + verDesc);

				for (SimplePathNodeEdge edge : path_overlap_graph.getOutEdges(vertex)) //get all edges of vertex->?
				{
					toVertex = path_overlap_graph.getDest(edge);

					p.println("\t" + vertex.getPathNodeID() + "->" + toVertex.getPathNodeID());

				}
			}


			p.println("}");
			
			p.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	

	}

	
	
	private static DirectedSparseGraph<Path, SimplePathNodeEdge> construct_path_overlap_graph(
			List<Path> path_list, HashMap<String, PathOverlap> pathMatches, HashSet<Integer> dispersed_repeat_nodes, String dot_file_prefix, String graphName, boolean createMiddleDotFiles) {
		
		// draw an edge between each pathNode B and the pathNode A to which B has a best-matching extension to the right.
		
		DirectedSparseGraph<Path, SimplePathNodeEdge> path_overlap_graph = new DirectedSparseGraph<Path, SimplePathNodeEdge>();
		
		for (Path p : path_list) {
			path_overlap_graph.addVertex(p);
		}
		
		// identify repeat nodes.
		HashSet<Integer> repeat_node_ids = new HashSet<Integer>();
		if (! dispersed_repeat_nodes.isEmpty()) {
			repeat_node_ids.addAll(dispersed_repeat_nodes);
		}
		
		for (Path path : path_list) {
			HashMap<Integer,Integer> repeat_nodes_and_counts = Path.getRepeatNodesAndCounts(path.get_vertex_list());
			for (Integer i : repeat_nodes_and_counts.keySet()) {
				repeat_node_ids.add(i);
			}
		}


		boolean store_best_extension_match_only = false;
		
		for (int i = 0; i < path_list.size(); i++) {
			
			int best_match = 0;
			int best_matching_path_idx = -1;
			
			List<Integer> best_precursor_j_indices = new ArrayList<Integer>();
			
			
			for (int j = 0; j < path_list.size(); j++) {
				
				if (i==j) {continue;}
				
				PathOverlap path_overlap = Path.pathB_extends_pathA_allowRepeats(path_list.get(i).get_vertex_list(),
											   								path_list.get(j).get_vertex_list(), 
											   								repeat_node_ids);
				
				
				int extension_matches = path_overlap.match_score;
				
				if (extension_matches <= 0) {
					continue;
				}
				
				// i extends j
				
				// got a match.
				String path_pair_token = get_path_compare_token(path_list.get(j), path_list.get(i));
				pathMatches.put(path_pair_token, path_overlap);

				debugMes("PathNode Overlap Detected: [overlap:  " + path_overlap.match_length + "] " 
						+ path_list.get(j) + " extended by " + path_list.get(i), 15);
				

				if (! store_best_extension_match_only)  {
					// add edge
					best_precursor_j_indices.add(j);

				}
				else {
					
					// examine for best extension

					if (extension_matches > best_match) {
						best_match = extension_matches;
						best_matching_path_idx = j;
						best_precursor_j_indices.clear();
						best_precursor_j_indices.add(j);


					}
					else if (extension_matches == best_match) {
						best_precursor_j_indices.add(j);
					}
				}
				
			}
			
		
		
			// add edges between overlapping and compatible paths:
			if (best_precursor_j_indices.size() > 0) {
				for (Integer precursor_index : best_precursor_j_indices) {

					String path_pair_token = get_path_compare_token(path_list.get(precursor_index), path_list.get(i));
					PathOverlap po = pathMatches.get(path_pair_token);

					
					debugMes("extension of: " + path_list.get(precursor_index) + " by " + path_list.get(i) 
							+  " has " + po.match_score + " terminal matches.", 15);

					// i extends j
					SimplePathNodeEdge spne = new SimplePathNodeEdge(po.match_score, 
							path_list.get(precursor_index).getPathNodeID(),
							path_list.get(i).getPathNodeID());
					
					path_overlap_graph.addEdge(spne, path_list.get(precursor_index), path_list.get(i));

				}
			}
			else {
				debugMes("path " + path_list.get(i) + " extends no path", 15);
			}
			
			
			
		}
		
		
		return(path_overlap_graph);
		
	}


	private static String get_path_compare_token(Path pathA, Path pathB) {
		
		String token = pathA.getPathNodeID() + ";" + pathB.getPathNodeID();
		
		return(token);
		
		
	}


	private static List<List<Integer>> remove_containments(
			List<List<Integer>> paths, HashMap<List<Integer>, 
			List<List<Integer>>> contained_path_to_containers) {
	
		// paths should already be sorted by descending length
		
		
		List<List<Integer>> noncontained_paths = new ArrayList<List<Integer>>();
		
		for (List<Integer> path : paths) {
			
			boolean contained = false;
			
			for (List<Integer> chosen_path : noncontained_paths) {
			
				if (Path.pathA_contains_pathB_allowRepeats(chosen_path, path)) {
					contained = true;
					
					// store containment info
					if (! contained_path_to_containers.containsKey(path)) {
						contained_path_to_containers.put(path, new ArrayList<List<Integer>>());
					}
					contained_path_to_containers.get(path).add(chosen_path);
				}
				
			}
			if (! contained) {
				noncontained_paths.add(path);
			}
			
		}
		
		return(noncontained_paths);
		
	}


	private static  DirectedSparseGraph<SeqVertex, SimpleEdge> construct_acyclic_graph(
			DirectedSparseGraph<SeqVertex, SimpleEdge> orig_graph, 
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(combinedReadHash, pairPaths, pairPathToReadSupport);
		
		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (PairPath pp : pairPaths) {
			
			paths.add(pp.getPath1());
			if (pp.hasSecondPath()) {
				paths.add(pp.getPath2());
			}
		}

		Collections.sort(paths, new Comparator<List<Integer>>() {

			public int compare(List<Integer> pathA, List<Integer> pathB) {

				if (pathA.size() < pathB.size()) {
					return(-1);
				}
				else if (pathA.size() > pathB.size()) {
					return(1);
				}
				else {
					return(0);
				}
			}



		});

		Collections.reverse(paths); // want descending by path l
		
		DirectedSparseGraph<SeqVertex, SimpleEdge> new_graph = 
				new DirectedSparseGraph<SeqVertex,SimpleEdge>();
		
		List<List<Integer>> cycle_inducing_paths = new ArrayList<List<Integer>>();
		HashSet<SimpleEdge> cycle_inducing_edges = new HashSet<SimpleEdge>();
		
		//DijkstraShortestPath<SeqVertex, SimpleEdge> dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(orig_graph);

		for (List<Integer> path : paths) {
			add_path_to_graph_disallow_cycles(orig_graph, new_graph, path, cycle_inducing_paths, cycle_inducing_edges);
		}
		
		debugMes("\n\nAll loop-inducing edges are: " + cycle_inducing_edges + "\n\ncontained in loop-inducing paths: " + cycle_inducing_paths, 10);
		
		
		return(new_graph);
	}


	private static void add_path_to_graph_disallow_cycles(
			DirectedSparseGraph<SeqVertex, SimpleEdge> orig_graph,
			DirectedSparseGraph<SeqVertex, SimpleEdge> new_graph, 
			List<Integer> path, 
			List<List<Integer>> cycle_inducing_paths,
			HashSet<SimpleEdge> cycle_inducing_edges) {

		debugMes("-adding path to new graph: " + path, 10);
		
		
		if (path.size() == 1) {
			SeqVertex v = SeqVertex.retrieveSeqVertexByID(path.get(0));
			new_graph.addVertex(v);
		}
		
		
		boolean cycle_inducing_path = false;
		for (int i = 1; i < path.size(); i++) {
			SeqVertex prev_node = SeqVertex.retrieveSeqVertexByID(path.get(i-1));
			SeqVertex next_node = SeqVertex.retrieveSeqVertexByID(path.get(i));
			
			boolean both_nodes_already_exist_in_graph = true;
			if (! new_graph.containsVertex(prev_node)) {
				new_graph.addVertex(prev_node);
				both_nodes_already_exist_in_graph = false;
			}
			if (! new_graph.containsVertex(next_node)) {
				new_graph.addVertex(next_node);
				both_nodes_already_exist_in_graph = false;
			}
			
			boolean add_edge = false;
			
			SimpleEdge se = orig_graph.findEdge(prev_node, next_node);
			
			if (both_nodes_already_exist_in_graph) {
			
				DijkstraShortestPath<SeqVertex, SimpleEdge> dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(new_graph);

				

				if (dp.getDistance(next_node, prev_node) != null) {
					// adding an edge between prev->next node would create a cycle!!!  

					cycle_inducing_path = true;
					
					cycle_inducing_edges.add(se);
					debugMes("\t** cycle-inducing edge found: " + prev_node + " to " + next_node, 10);

				}
				else {
					// add edge to graph.
					add_edge = true;
				}
			}
			else {
				add_edge = true;
			}
			
			if (add_edge) {
				new_graph.addEdge(se, prev_node, next_node);
			}
		}
		
		if (cycle_inducing_path) {
			cycle_inducing_paths.add(path);
			debugMes("\t$$ cycle inducing path: " + path, 15);
		}
		
		
	}


	private static boolean graph_contains_loops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		
		DijkstraShortestPath<SeqVertex, SimpleEdge> dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(graph);

		///////////////////////////////////////////////////////////////////////////////////////
		//  ===== Capture all remaining loops found. ========== 
		//
		// These should be only those repeats that aren't evident in the individual read paths,
		// since the read-evident repeats were unrolled earlier.
		/////////////////////////////////////////////////////////////////////////////////////////

		Set<Set<SimpleEdge>> curLoops = new HashSet<Set<SimpleEdge>>();

		// find all loops in the graph by seeing if, given edge v->v2, there is a path from v2 back to v
		for (SeqVertex v : graph.getVertices())
		{
			for (SeqVertex v2 : graph.getSuccessors(v))
			{
				if (dp.getDistance(v2, v)!=null) // there is a connection between v->v2->... ->v
				{
					//path has all edges from v to itself thru v2
					List<SimpleEdge> loopPath = dp.getPath(v2, v);

					// v2 is successor of v, so let's just add the v->v2 edge too, complete the full loop.
					loopPath.add(0, graph.findEdge(v, v2)); 

					// capture the path IDs for debugMes reporting below.
					List<Integer> pathIDs = new ArrayList<Integer>();
					for (SimpleEdge e : loopPath)
						pathIDs.add(graph.getDest(e).getID());

					// Collect the loop edge set.
					Set<SimpleEdge> loopPath_set = new HashSet<SimpleEdge>(loopPath);

					if (!curLoops.contains(loopPath_set))
					{
						curLoops.add(loopPath_set);
						debugMes("adding the loop path "+pathIDs+" to the curLoops",12);
					}else
					{
						debugMes("not adding the loop path "+pathIDs+" to the curLoops",12);
					}
				}
			}
		}

		if (curLoops.isEmpty()) {
			return(false);
		}
		else {
			return(true);
		}


	}


	private static HashMap<Integer, HashMap<PairPath, Integer>> reassign_read_paths_according_to_longer_path_compatibility(
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		
		debugMes("\nSECTION\n========= Ressigning Repeat-containing Read Paths Based On Longer Path Compatibility ==========\n\n",5);
		
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(combinedReadHash, pairPaths, pairPathToReadSupport);
		
		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (PairPath pp : pairPaths) {
			
			paths.add(pp.getPath1());
			if (pp.hasSecondPath()) {
				paths.add(pp.getPath2());
			}
		}

		Collections.sort(paths, new Comparator<List<Integer>>() {

			public int compare(List<Integer> pathA, List<Integer> pathB) {

				if (pathA.size() < pathB.size()) {
					return(-1);
				}
				else if (pathA.size() > pathB.size()) {
					return(1);
				}
				else {
					return(0);
				}
			}



		});

		Collections.reverse(paths); // want descending by path l
		
		// convert to list of PathWithOrig
		List<PathWithOrig> path_with_orig_list = new ArrayList<PathWithOrig>();
		
		for (List<Integer> path : paths) {
			path_with_orig_list.add(new PathWithOrig(path));
			
		}
		
		
		// iterate through pair paths and see if they require reassignment.
		HashMap<PairPath,Integer> updated_pairpath_hmap = new HashMap<PairPath,Integer>();
		
		for (PairPath pp : pairPaths) {
			
			Integer read_support = pairPathToReadSupport.get(pp);
			
			PairPathWOrig ppwo = new PairPathWOrig(pp);
			
			boolean restructured_flag = false;
			
			for (PathWithOrig template_pwo : path_with_orig_list) {
				
				if (template_pwo.size() < ppwo.size()){
					break;
				}
				
				PairPathWOrig ppwo_restructured = ppwo.restructure_according_to_repeat_path(template_pwo);
				if (ppwo_restructured != null) {
					updated_pairpath_hmap.put(ppwo_restructured.getPairPath(), read_support);
					
					if (! ppwo_restructured.equals(ppwo)) {
						debugMes("PPWO restructured from: " + ppwo + " to " + ppwo_restructured, 15);
					}
					
					restructured_flag = true;
					break;
				}
			}
			if (! restructured_flag) {
				// stick with the original one.
				
				if (true) {
					throw new RuntimeException("error, not restructured: " + path_with_orig_list + " and target: " + ppwo);
				}
				updated_pairpath_hmap.put(pp, read_support);
				
			}
		}
		
		HashMap<Integer, HashMap<PairPath, Integer>> new_combined_read_hash = construct_combinedReadhHash_from_PairPath_list(updated_pairpath_hmap);
		
		return(new_combined_read_hash);
		
		
	}


	private static HashMap<Integer, HashMap<PairPath, Integer>> construct_combinedReadhHash_from_PairPath_list(
			HashMap<PairPath, Integer> pairpath_hmap) {
		
		HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash = new HashMap<Integer, HashMap<PairPath, Integer>>();
		
		for (PairPath pp : pairpath_hmap.keySet()) {
			
			Integer read_support = pairpath_hmap.get(pp);
			
			Integer first_id = pp.getFirstID();
			
			HashMap<PairPath,Integer> pp_map;
			
			if (combinedReadHash.containsKey(first_id)) {
				pp_map = combinedReadHash.get(first_id);
			}
			else {
				pp_map = new HashMap<PairPath,Integer>();
				combinedReadHash.put(first_id, pp_map);	
			}
			
			pp_map.put(pp, read_support);
			
		}
		
		return(combinedReadHash);
	}


	private static void examine_out_of_order_depth_in_read_paths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash) {
	
		
		debugMes("\n\nSECTION\n====  examining node depths of read paths in DAG ======\n\n", 5);
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(componentReadHash, pairPaths, pairPathToReadSupport);
		
		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (PairPath pp : pairPaths) {
			paths.add(pp.getPath1());
			if (pp.hasSecondPath()) {
				paths.add(pp.getPath2());
			}
		}

		Collections.sort(paths, new Comparator<List<Integer>>() {

			public int compare(List<Integer> pathA, List<Integer> pathB) {

				if (pathA.size() < pathB.size()) {
					return(-1);
				}
				else if (pathA.size() > pathB.size()) {
					return(1);
				}
				else {
					return(0);
				}
			}



		});

		Collections.reverse(paths); // want descending by path length

		
		int count_reads_ok = 0;
		int count_reads_conflict = 0;
		
		for (List<Integer> path : paths) {

			List<List<Integer>> read_parts = new ArrayList<List<Integer>>();

			List<Integer> part = new ArrayList<Integer>();

			List<List<Integer>> node_depths_tracker = new ArrayList<List<Integer>>();



			ArrayList<Integer> node_depths_list = new ArrayList<Integer>();

			Iterator<Integer> it = path.iterator();
			int prev_depth = -1;

			HashMap<Integer,Boolean> node_visitor = new HashMap<Integer,Boolean>();
			while (it.hasNext()) {
				Integer node_id = it.next();

				SeqVertex v = getSeqVertex(graph, node_id);
				if (v._node_depth < prev_depth || node_visitor.containsKey(node_id)) {
					// problem...
					// fracture here.
					if (! part.isEmpty()) {
						read_parts.add(part);
						part = new ArrayList<Integer>();

					}
				}
				node_depths_list.add(v._node_depth);
				prev_depth = v._node_depth;
				node_visitor.put(node_id, true);
				part.add(node_id);
			}
			if (! part.isEmpty()) {
				read_parts.add(part);
				part = new ArrayList<Integer>();
			}

			node_depths_tracker.add(node_depths_list);
		
			if (read_parts.size() > 1) {
			
				debugMes("DAG-conflicting path: " + path + " with node_depths: " + node_depths_tracker  + " into " + read_parts.size() + ": " + read_parts, 10);
				
				count_reads_conflict++;
			}
			else {
				count_reads_ok++;
			}

		}
		
		debugMes("\n\nNum reads ok: " + count_reads_ok + "\nNum reads conflicted depths: " + count_reads_conflict, 10);
		
	}


	private static void examine_compatible_paths_debugging_only(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash) {
		
		// this method is only for debugging purposes.	
		
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(componentReadHash, pairPaths, pairPathToReadSupport);
		
		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (PairPath pp : pairPaths) {
			paths.add(pp.getPath1());
			if (pp.hasSecondPath()) {
				paths.add(pp.getPath2());
			}
		}

		Collections.sort(paths, new Comparator<List<Integer>>() {

			public int compare(List<Integer> pathA, List<Integer> pathB) {

				if (pathA.size() < pathB.size()) {
					return(-1);
				}
				else if (pathA.size() > pathB.size()) {
					return(1);
				}
				else {
					return(0);
				}
			}



		});

		Collections.reverse(paths); // want descending by path length


		List<List<Integer>> longest_paths = new ArrayList<List<Integer>>();
		
		////////////////////////////////
		int NUM_TOP_LONGEST_PATHS = 2; // change to the number of longest paths to compare to.
		///////////////////////////////
		
		
		int counter = 0;
		
		for (List<Integer> path : paths) {
			counter++;
			if (counter <= NUM_TOP_LONGEST_PATHS) {
				longest_paths.add(path);
			}
			else {
				boolean compatible = false;

				for (List<Integer> longer_path : longest_paths) {

					if (PairPath.individual_paths_are_compatible(path, longer_path)) {
						compatible = true;
						break;
					}
				}
				
				if (! compatible) {
					debugMes("LongPathIncompat: " + path, 10);
				
				}
			}
			
		}
		
		
	}


	private static HashMap<List<Integer>, Pair<Integer>> reconstruct_paths_from_collapsed_reads(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash) {
HashMap<List<Integer>, Pair<Integer>> transcripts = new HashMap<List<Integer>,Pair<Integer>>();
		
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(componentReadHash, pairPaths, pairPathToReadSupport);
		
		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (PairPath pp : pairPaths) {
			paths.add(pp.getPath1());
			if (pp.hasSecondPath()) {
				paths.add(pp.getPath2());
			}
		}
		
		Collections.sort(paths, new Comparator<List<Integer>>() {
			
			public int compare(List<Integer> pathA, List<Integer> pathB) {
				
				if (pathA.size() < pathB.size()) {
					return(-1);
				}
				else if (pathA.size() > pathB.size()) {
					return(1);
				}
				else {
					return(0);
				}
			}
			
			
			
		});
		
		Collections.reverse(paths); // want descending by path length
		
		
		List<List<Integer>> collapsed_paths = Path.collapse_compatible_paths_to_min_set(paths);
		
		for (List<Integer> path : collapsed_paths) {
			
			// only doing pp path1 since testing for cufflinks here in unpaired mode.
			transcripts.put(path, new Pair(new Integer(1), new Integer(1)));
		
		}
		
		return(transcripts);
	}


	private static HashMap<List<Integer>, Pair<Integer>> reconstruct_paths_from_reads(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash) {

		
		HashMap<List<Integer>, Pair<Integer>> transcripts = new HashMap<List<Integer>,Pair<Integer>>();
		
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(componentReadHash, pairPaths, pairPathToReadSupport);
		
		for (PairPath pp : pairPaths) {
			
			// only doing pp path1 since testing for cufflinks here in unpaired mode.
			transcripts.put(pp.getPath1(), new Pair(new Integer(1), new Integer(1)));
		
		}
		
		return(transcripts);
		
	}


	private static HashMap<List<Integer>, Pair<Integer>> remove_lower_ranked_paths_without_unique_read_content(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, HashMap<List<Integer>, Pair<Integer>> finalPaths_all,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> finalPathsToContainedReads) {
		
		
		
		
		HashMap<PairPath,Boolean> pp_used = new HashMap<PairPath,Boolean>();
		
		
		HashSet<List<Integer>> all_paths = new HashSet<List<Integer>>();
		
		// init the list of candidate paths
		List<List<Integer>> paths = new ArrayList<List<Integer>>();
		for (List<Integer> path : finalPaths_all.keySet()) {
			paths.add(path);
			all_paths.add(path);
		}
		
		HashMap<List<Integer>,Integer> seqLengthMap = paths_to_seq_lengths(graph, all_paths);
		
		HashSet<List<Integer>> priority_paths_with_unique_read_content = new HashSet<List<Integer>>();
		
		

		int round = 0;
		
		while (! paths.isEmpty()) {
			
			round++;
			
			UniquePathContentComparator unique_path_content_comparator = new UniquePathContentComparator(paths, 
																			pp_used, finalPathsToContainedReads,
																			seqLengthMap);
			
			paths = unique_path_content_comparator.remove_paths_without_unique_read_content(paths);
			
			if (! paths.isEmpty()) {
				

				Collections.sort(paths, unique_path_content_comparator);
				Collections.reverse(paths);

				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
					debugMes("Round[" + round + "] Paths Prioritized by Unique Read Content", 15);
					for (List<Integer> path : paths) {
						debugMes("Round[" + round + "] Unique=" + unique_path_content_comparator.unique_count(path) + ", path: " + path, 15);
					}
					
					debugMes("Round[ " + round + "] SELECTING Unique=" + unique_path_content_comparator.unique_count(paths.get(0)) 
							+ ", path: "  + paths.get(0) + "\n\n", 15);
				}
				
				
				List<Integer> p = paths.remove(0);
				
				priority_paths_with_unique_read_content.add(p);
			
				// mark contained paths as seen
				for (PairPath pp : finalPathsToContainedReads.get(p).keySet()) {
					pp_used.put(pp, true);
				}
			}
			
		}
		
		
		HashMap<List<Integer>, Pair<Integer>> paths_to_keep = new HashMap<List<Integer>, Pair<Integer>>();
		
		for (List<Integer> path : all_paths) {
			paths_to_keep.put(path, finalPaths_all.get(path));
			debugMes("PathRankingFilter retaining: " + path, 15);
		}
		
		return(paths_to_keep);
		
		
	}


	private static HashMap<List<Integer>, Integer> paths_to_seq_lengths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashSet<List<Integer>> all_paths) {
		
		
		HashMap<List<Integer>, Integer> seqLengthMap = new HashMap<List<Integer>,Integer>();


		for (List<Integer> path : all_paths)
		{
			String seq = getPathSeq(graph,path);

			seqLengthMap.put(path, seq.length());



		}

		return(seqLengthMap);
	}
		
	
	
	private static void removeShortOrphanNodes(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			int min_seq_length) {
		
		List<SeqVertex> vertices_to_remove = new ArrayList<SeqVertex>();

		for (SeqVertex v : graph.getVertices()) {

			if (graph.getPredecessorCount(v) == 0
				&&
				graph.getSuccessorCount(v) == 0
				&&
				v.getName().length() < min_seq_length) 
			{
					vertices_to_remove.add(v);
			}
		}
		
		for (SeqVertex v : vertices_to_remove) {
			debugMes("Removing short seq orphaned vertex: " + v + " from graph. Seq too short to generate a contig of min length.", 12);
			graph.removeVertex(v);
			
		}
		
	}


	private static HashMap<Integer, HashMap<PairPath, Integer>> reorganizeReadPairings(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis) {
		
		
		describeNodes(graph);
		
		
		
		HashMap<Integer,HashMap<PairPath,Integer>>   newCombinedReadHash = new HashMap<Integer,HashMap<PairPath,Integer>>  ();

		for (HashMap<PairPath, Integer> pairs_n_counts : combinedReadHash.values()) {
			
			for (PairPath pp : pairs_n_counts.keySet()) {
				
				int read_support = pairs_n_counts.get(pp);
				
				if (pp.hasSecondPath()) {

					PairPath  combinedPath = combinePaths(graph, pp.getPath1(), pp.getPath2(), dijkstraDis);

					if (! combinedPath.isEmpty()) {
						storePairPathByFirstVertex(combinedPath, newCombinedReadHash, read_support);
						debugMes("OK pp update to new DAG: " + pp + " => " + combinedPath, 15);
					}
					else {
						// store the read path separately
						PairPath pp1 = new PairPath(pp.getPath1());
						storePairPathByFirstVertex(pp1, newCombinedReadHash, read_support);

						PairPath pp2 = new PairPath(pp.getPath2());
						storePairPathByFirstVertex(pp2, newCombinedReadHash, read_support);
						debugMes("Warning... pp: " + pp + " needed to be split into: " + pp1 + " and " + pp2, 15);
					}
				}
				else {
					storePairPathByFirstVertex(pp, newCombinedReadHash, read_support);
				}

			}
		}
		
		
		return(newCombinedReadHash);
		
	}


	private static void storePairPathByFirstVertex(PairPath pp,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash, int read_support) {
		
		Integer firstV = pp.getFirstID();

		if (!combinedReadHash.containsKey(firstV))
			combinedReadHash.put(firstV, new HashMap<PairPath,Integer>()); //init

		if (!combinedReadHash.get(firstV).containsKey(pp))
			combinedReadHash.get(firstV).put(pp, 0); //add pairpath

		Integer counts = combinedReadHash.get(firstV).get(pp);
		combinedReadHash.get(firstV).put(pp, counts + read_support); // increment counts for pairpath
		debugMes("we have "+ combinedReadHash.get(firstV).get(pp)+" reads supporting the path: " + pp,18);
		
		
	}
	

	private static int unroll_remaining_terminal_self_loops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		
		int num_self_loops_unrolled = 0;
		
		List<SeqVertex> all_vertices = new ArrayList<SeqVertex>(graph.getVertices());
		
		HashMap<Integer,Boolean> unrolled_terminal_vertices = new HashMap<Integer,Boolean>();
		
		for (SeqVertex v : all_vertices) {
			
			debugMes("Examining node: " + v.getShortSeqWconnectingIDs(graph), 12);
			
			Collection<SeqVertex> successors = graph.getSuccessors(v);
			Collection<SeqVertex> preds = graph.getPredecessors(v);
			
			if (successors.size() == 1 && preds.size() == 1 
					&& successors.containsAll(preds)) {
			
				// a terminal repeat node. 
				
				
				
				SeqVertex repeat_vertex = successors.iterator().next();
				
				if (repeat_vertex.getOrigButterflyID() != v.getOrigButterflyID()) {
					debugMes("not a terminal self loop, skipping..." 
							+ v, 12);
					continue;
				}
				
				debugMes("Removing terminal self loop at vertex: " + v, 12);
			
				SeqVertex new_v = new SeqVertex(getNextID(), repeat_vertex); // this constructor sets orig_id so it's the same.

				SimpleEdge loop_edge = graph.findEdge(v, repeat_vertex);
				
				graph.addVertex(new_v);
				
				SimpleEdge new_edge = new SimpleEdge(loop_edge.getWeight(), v.getID(), new_v.getID());
				graph.addEdge(new_edge, v, new_v);
				new_edge.increment_repeat_unroll_weight(2);
				
				num_self_loops_unrolled++;
				
				unrolled_terminal_vertices.put(v.getID(), true);
				unrolled_terminal_vertices.put(repeat_vertex.getID(), true);  // flag for possible read reassignment
				
				// copy over any outgoing non-self edges from the repeat vertex
				for (SeqVertex succ : graph.getSuccessors(repeat_vertex)) {
					if (succ.getOrigButterflyID() != new_v.getOrigButterflyID()) {
						
						SimpleEdge se = graph.findEdge(repeat_vertex, succ);
						SimpleEdge new_se = new SimpleEdge(se.getWeight(), new_v.getID(), succ.getID());
						graph.addEdge(new_se, new_v, succ);
						new_se.increment_repeat_unroll_weight(2);
						unrolled_terminal_vertices.put(succ.getID(), true); // flag for possible reassignment
						
					}
				}
				
				
			}
		
			
		}
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
			try {
				writeDotFile(graph, "__terminal_loop_unroll-preReadReassignment.dot", "ladeda", false);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		
		
		int num_paths_redefined = reassign_repeat_nodes_in_reads(graph, 
				combinedReadHash, unrolled_terminal_vertices, 
				null, false);
		
		debugMes("num paths with terminal self loop vertices and paths redefined: " + num_paths_redefined, 12);
		
		/*
		if (num_paths_redefined > 0)
			redefine_all_graph_edges(graph, combinedReadHash);
		*/
		
		
		return(num_self_loops_unrolled);
	}


	private static void describeVertices(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		debugMes("## Node descriptions:", 10);
		
		List<SeqVertex> vertices = new ArrayList<SeqVertex>(graph.getVertices());
		

		Collections.sort(vertices, SeqVertexIDorderComparator);
		
		for (SeqVertex v : vertices) {
			debugMes(v.getShortSeqWconnectingIDs(graph), 10);
		}
		
	}


	private static void describeNodes(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		for (SeqVertex v : graph.getVertices()) {
			debugMes("NODE_DESCR: " + v.getShortSeqWconnectingIDs(graph), 5);
		}
		
		
	}


	private static boolean unroll_loops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash, My_DFS dfs) {
	
		debugMes("\n\nUNROLLING LOOPS IN READS\n\n", 5);
		
		
		HashMap<Integer,Boolean> all_repeat_related_nodes =  new HashMap<Integer,Boolean>();
		
		// get the repeat nodes assigned to the pp containing it the greatest number of times.
		HashMap<Integer, PairPath> repeat_node_id_to_longest_path = find_repeat_containing_pairpaths_ignoreLastNode(combinedReadHash);
		
		if (repeat_node_id_to_longest_path.isEmpty()) {
			debugMes("\t** no repeats detected in the reads. No repeat unrolling needed here.", 10);
			return(false);
		}
		
		
		
		int unroll_loop_counter = 0;
		
		while (repeat_node_id_to_longest_path.size() > 0) {
			
			unroll_loop_counter++;
			
			debugMes("\n\n## Unrolling loops, round: " + unroll_loop_counter, 10);
			
			
			if (repeat_node_id_to_longest_path.size() > 0) {
				debugMes("\n\nFound : " + repeat_node_id_to_longest_path.size() + " repeat nodes.", 10);
				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 12) {
					for (Integer node_id : repeat_node_id_to_longest_path.keySet()) {
						PairPath pp = repeat_node_id_to_longest_path.get(node_id);
						System.err.println("Repeat Node: " + node_id + " found in longest pp: " + pp);
					}
					
				}
				
			}
			
			
			
			// unroll each repeat unit
			HashMap<Integer,Boolean> restructured_nodes = new HashMap<Integer,Boolean>();
		
			
			// prioritize the paths according to the number of different repeat nodes.
			HashMap<PairPath,Float> pp_uniq_repeat_count_hmap = new HashMap<PairPath,Float>();
			for (Integer repeat_node : repeat_node_id_to_longest_path.keySet()) {
				
				all_repeat_related_nodes.put(repeat_node, true);
				
				PairPath pp = repeat_node_id_to_longest_path.get(repeat_node);
				
				float repeat_count_sum = pp.getRepeatNodesAndCountSum();
				// also prioritize by path length
				int max_path_length = pp.getMaxPathLength();
				
				repeat_count_sum += (float)max_path_length/1000.0; // simple way to augment score based on length of longest path.
				
				if (! pp_uniq_repeat_count_hmap.containsKey(pp)) {
					pp_uniq_repeat_count_hmap.put(pp, repeat_count_sum); // all repeat nodes counted, do once per pp
				}
				
			}
			
			final HashMap<PairPath,Float> pp_uniq_repeat_count_hmap_copy = new HashMap<PairPath,Float>(pp_uniq_repeat_count_hmap);
			
			List<PairPath> pp_list_ordered_by_repeat_counts = new ArrayList<PairPath>(pp_uniq_repeat_count_hmap.keySet());
			Collections.sort(pp_list_ordered_by_repeat_counts, new Comparator<PairPath>() {

				public int compare (PairPath pp_A, PairPath pp_B) {

					float rpt_node_count_A = pp_uniq_repeat_count_hmap_copy.get(pp_A);
					float rpt_node_count_B = pp_uniq_repeat_count_hmap_copy.get(pp_B);

					if (rpt_node_count_A > rpt_node_count_B)
						return -1;
					else if (rpt_node_count_A < rpt_node_count_B)
						return 1;
					else
						return 0;
					
				}

			});


			for (PairPath pp : pp_list_ordered_by_repeat_counts) {

				debugMes("Unrolling repeats in pp: " + pp, 12);

				for (List<Integer> path : pp.get_paths()) {

					if (path.size() == 0) // empty path in a pp
						continue;


					// if the path contains any previously restructured node, skip it.
					if (Path.contains_any_node_id(path, restructured_nodes.keySet())) {
						debugMes("\t-postponing unroll since contains restructured node.\n", 12);
						continue;
					}

					debugMes("Unrolling repeats in pp: " + pp + " with repeat nodes: " + pp.getRepeatNodesAndCounts(), 12);

					// find the repeat nodes and restructure the graph.
					List<SeqVertex> path_vertices = new ArrayList<SeqVertex>();
					HashMap<Integer,Boolean> seen_repeat_node = new HashMap<Integer,Boolean>();

					boolean restructured_flag = false;

					int node_pos = 0;
					for (Integer node_id : path) {

						node_pos++;

						SeqVertex v = getSeqVertex(graph, node_id);

						if (node_pos != path.size())  {
							// don't do the very last node if it's a repeat... leave for another pp to unravel it.

							if (seen_repeat_node.containsKey(node_id)) {

								// ok - tackling this repeat node.  

								//////////////////////////////
								// make a copy of the vertex
								//////////////////////////////

								v = new SeqVertex(getNextID(), v); // this constructor sets orig_id so it's the same.

								debugMes("\tcopying node: " + node_id + " to " + v.getID(), 12);

								all_repeat_related_nodes.put(v.getID(), true);

								graph.addVertex(v);
								restructured_nodes.put(node_id, true);

								restructured_flag = true;

							}
							else {


								if (repeat_node_id_to_longest_path.containsKey(node_id)) {
									seen_repeat_node.put(node_id, true);
								}
							}



						} // endif last node of path

						path_vertices.add(v);

					} //end of for node_id : path


					if (restructured_flag) {

						// add edges between neighboring nodes
						for (int i = 0; i < path_vertices.size()-1; i++) {
							SeqVertex before_node = path_vertices.get(i);
							SeqVertex after_node = path_vertices.get(i+1);

							SimpleEdge edge = graph.findEdge(before_node, after_node);
							if (edge == null) {
								// add it.
								// get the original edge and copy the weight over.
								SimpleEdge orig_edge = graph.findEdge(getSeqVertex(graph, before_node.getOrigButterflyID()), getSeqVertex(graph, after_node.getOrigButterflyID()));
								double oldW = 1;
								if (orig_edge != null) 
									oldW = orig_edge.getWeight();

								edge = new SimpleEdge(oldW, before_node.getID(), after_node.getID());

								graph.addEdge(edge, before_node, after_node);

							}
							edge.increment_repeat_unroll_weight(unroll_loop_counter); // was 1


						}

						// describe the new vertex list:
						if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
							debugMes("# Restructured path described:", 20);
							for (SeqVertex v : path_vertices) {
								debugMes(v.getShortSeqWconnectingIDs(graph), 20);
							}


						}




						// Verify that this path can be properly reassigned in the graph.

						debugMes("\nVerifying that restructured path: " + path + " is rethreaded through the graph with fewer repeat units.", 12);
						List<Integer> updated_path = reassign_restructured_path_in_graph(graph, path);

						String orig_path_seq = getPathSeq(graph, path);
						String new_path_seq = getPathSeq(graph, updated_path);
						
						if (! orig_path_seq.equals(new_path_seq)) {
							throw new RuntimeException("Error, updated path seq != orig path seq:\n>Orig\n" + orig_path_seq + "\n>New\n" + new_path_seq);
						}
						else {
							debugMes("* old and new path seqs are identical. validated. " + new_path_seq, 15);
						}


						if (updated_path == null 
								|| 

								(   path.size() > 1 &&	
										(updated_path.equals(path)
												||
												Path.countNumNodesNotUnique(path) <= Path.countNumNodesNotUnique(updated_path) 
												)

										)

								) {

							throw new RuntimeException("Repeat-unrolled path: " + path 
									+ " " + Path.getRepeatNodesAndCounts(path)
									+ " was not properly restructured: "
									+ updated_path + " " + Path.getRepeatNodesAndCounts(updated_path));

						}
						else {
							debugMes("\tVerification OK: path:" + path + " " + Path.getRepeatNodesAndCounts(path) + 
									" => " +  updated_path + " " + Path.getRepeatNodesAndCounts(updated_path) + "\n", 12); 
						}

					} // end if restructured_flag
					
					
				} // end for path
				
				
			} // end for pp
			
			
			
			
			
			// reassign repeat nodes to their new nodes in the graph
			debugMes("\n\n## Post-unroll round: " + unroll_loop_counter + ", reassigning_repeat_nodes_in_reads\n", 10);

			// restrict unrolling to just those containing as of yet unrolled repeats
			int num_paths_updated = reassign_repeat_nodes_in_reads(graph, combinedReadHash, 
					restructured_nodes, null, true);

			if (num_paths_updated == 0)
				throw new RuntimeException("Error, no paths were updated after this round of repeat unrolling");
			
			debugMes("\n\nNumber of paths refined: " + num_paths_updated, 10);
			
			//redefine_all_graph_edges(graph, combinedReadHash);  // prune out the now unsupported edges post reassignment.
			

			dfs.runDFS2(); // reassign depths given new nodes added.

			if (true) {
				String filename = FILE + ".repeatUnroll_" + unroll_loop_counter + ".dot";
				PrintStream p;
				try {
					p = new PrintStream(new FileOutputStream(filename));
					writeDotFile(graph,p,"repeatUnroll_" + unroll_loop_counter);
					p.close();
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
			
			
			// look for remaining repeats
			repeat_node_id_to_longest_path = find_repeat_containing_pairpaths_ignoreLastNode(combinedReadHash);
			
			
		
		} // end while repeats
		
		// all repeats are now unrolled.
		
		
		
		
		// further refine paths in case better scoring path for a read is to be found.
		
		int num_paths_redefined = 1;

		HashMap<PairPath,Boolean> pp_remains_unchanged_skip_list = new HashMap<PairPath,Boolean>(); // now we use it. Graph isn't changing, so wont' expect paths to change unless they have a better placement given the unrolled edge weights.

		int refinement_round = 0;
		
		while (num_paths_redefined > 0) {

			refinement_round++;

			// reassign paths for all reads containing restructured nodes.
			num_paths_redefined = reassign_repeat_nodes_in_reads(graph, 
					combinedReadHash, all_repeat_related_nodes, 
					pp_remains_unchanged_skip_list, false);

			debugMes("unroll_loops::PATH_REFINEMENT_ROUND: " + refinement_round + " NUMBER_PATHS_REDEFINED: " + num_paths_redefined, 10);

			/*

			if (num_paths_redefined > 0 
					 || refinement_round == 1
			  )
			{ // be sure to do this at least once!!!
				redefine_all_graph_edges(graph, combinedReadHash);
			}
		
			*/
		
		}
		
		

		return(true);
	}


	
	private static boolean unroll_remaining_terminal_loops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
	
		// should be targeting only those repeat nodes that show up at final positions in a read.
		
		debugMes("\n\nUNROLLING REMAINING LOOPS IN READS\n\n", 5);
		
		
		HashMap<Integer,Boolean> all_repeat_related_nodes =  new HashMap<Integer,Boolean>();
		
		// get the repeat nodes assigned to the pp containing it the greatest number of times.
		HashMap<Integer, PairPath> repeat_node_id_to_longest_path = find_repeat_containing_pairpaths(combinedReadHash);
		
		if (repeat_node_id_to_longest_path.isEmpty()) {
			debugMes("\t** no repeats detected in the reads. No repeat unrolling needed here.", 10);
			return(false);
		}
		
		
		
		int unroll_loop_counter = 0;
		
		while (repeat_node_id_to_longest_path.size() > 0) {
			
			unroll_loop_counter++;
			
			debugMes("\n\n## Unrolling remaining terminal loops, round: " + unroll_loop_counter, 10);
			
			if (repeat_node_id_to_longest_path.size() > 0) {
				debugMes("\n\nFound : " + repeat_node_id_to_longest_path.size() + " repeat nodes.", 10);
				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 12) {
					for (Integer node_id : repeat_node_id_to_longest_path.keySet()) {
						PairPath pp = repeat_node_id_to_longest_path.get(node_id);
						System.err.println("Repeat Node: " + node_id + " found in longest pp: " + pp);
					}
					
				}
				
			}
			
			
			// unroll each repeat unit
			HashMap<Integer,Boolean> restructured_nodes = new HashMap<Integer,Boolean>();
		
			
			// prioritize the paths according to the number of different repeat nodes.
			HashMap<PairPath,Float> pp_uniq_repeat_count_hmap = new HashMap<PairPath,Float>();
			for (Integer repeat_node : repeat_node_id_to_longest_path.keySet()) {
				
				all_repeat_related_nodes.put(repeat_node, true);
				
				PairPath pp = repeat_node_id_to_longest_path.get(repeat_node);
				
				float repeat_count_sum = pp.getRepeatNodesAndCountSum();
				// also prioritize by path length
				int max_path_length = pp.getMaxPathLength();
				
				repeat_count_sum += (float)max_path_length/1000.0; // simple way to augment score based on length of longest path.
				
				if (! pp_uniq_repeat_count_hmap.containsKey(pp)) {
					pp_uniq_repeat_count_hmap.put(pp, repeat_count_sum); // all repeat nodes counted, do once per pp
				}
				
			}
			
			final HashMap<PairPath,Float> pp_uniq_repeat_count_hmap_copy = new HashMap<PairPath,Float>(pp_uniq_repeat_count_hmap);
			
			List<PairPath> pp_list_ordered_by_repeat_counts = new ArrayList<PairPath>(pp_uniq_repeat_count_hmap.keySet());
			Collections.sort(pp_list_ordered_by_repeat_counts, new Comparator<PairPath>() {

				public int compare (PairPath pp_A, PairPath pp_B) {

					float rpt_node_count_A = pp_uniq_repeat_count_hmap_copy.get(pp_A);
					float rpt_node_count_B = pp_uniq_repeat_count_hmap_copy.get(pp_B);

					if (rpt_node_count_A > rpt_node_count_B)
						return -1;
					else if (rpt_node_count_A < rpt_node_count_B)
						return 1;
					else
						return 0;
					
				}

			});


			for (PairPath pp : pp_list_ordered_by_repeat_counts) {

				debugMes("Unrolling repeats in pp: " + pp, 12);

				for (List<Integer> path : pp.get_paths()) {

					if (path.size() == 0) // empty path in a pp
						continue;


					// if the path contains any previously restructured node, skip it.
					if (Path.contains_any_node_id(path, restructured_nodes.keySet())) {
						debugMes("\t-postponing unroll since contains restructured node.\n", 12);
						continue;
					}

					debugMes("Unrolling repeats in pp: " + pp + " with repeat nodes: " + pp.getRepeatNodesAndCounts(), 12);

					// find the repeat nodes and restructure the graph.
					List<SeqVertex> path_vertices = new ArrayList<SeqVertex>();
					HashMap<Integer,Boolean> seen_repeat_node = new HashMap<Integer,Boolean>();

					boolean restructured_flag = false;

					int node_pos = 0;
					for (Integer node_id : path) {

						node_pos++;

						SeqVertex v = getSeqVertex(graph, node_id);


						if (seen_repeat_node.containsKey(node_id)) {

							// ok - tackling this repeat node.  

							
							// this had better be the very last node of the path, or else our logic is wrong here!!
							if (node_pos != path.size())
								throw new RuntimeException("Error, remaining repeat is not at last node of path: " + path);
							
							all_repeat_related_nodes.put(node_id, true); // important, track the node ID getting restructured.
							
							//////////////////////////////
							// make a copy of the vertex
							//////////////////////////////

							v = new SeqVertex(getNextID(), v); // this constructor sets orig_id so it's the same.

							debugMes("\tcopying node: " + node_id + " to " + v.getID(), 12);

							

							graph.addVertex(v);
							restructured_nodes.put(node_id, true);

							restructured_flag = true;

						}
						else {


							if (repeat_node_id_to_longest_path.containsKey(node_id)) {
								seen_repeat_node.put(node_id, true);
							}
						}
						
						
						path_vertices.add(v);

					} //end of for node_id : path


					if (restructured_flag) {

						// add edges between neighboring nodes
						for (int i = 0; i < path_vertices.size()-1; i++) {
							SeqVertex before_node = path_vertices.get(i);
							SeqVertex after_node = path_vertices.get(i+1);

							SimpleEdge edge = graph.findEdge(before_node, after_node);
							if (edge == null) {
								
								// again, this had better be an edge to the last node in the path
								if (i + 2 != path.size())
									throw new RuntimeException("Error, trying to add new edge between " 
											+ before_node + " and " + after_node + " and not at end of path.");
								
								// add it.
								// get the original edge and copy the weight over.
								SimpleEdge orig_edge = graph.findEdge(getSeqVertex(graph, before_node.getOrigButterflyID()), getSeqVertex(graph, after_node.getOrigButterflyID()));
								double oldW = 1;
								if (orig_edge != null) 
									oldW = orig_edge.getWeight();

								edge = new SimpleEdge(oldW, before_node.getID(), after_node.getID());

								graph.addEdge(edge, before_node, after_node);

								edge.increment_repeat_unroll_weight(2); // further encourage reads to take this path where possible.
							}
							
							

						}

						// describe the new vertex list:
						if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
							debugMes("# Restructured path described:", 20);
							for (SeqVertex v : path_vertices) {
								debugMes(v.getShortSeqWconnectingIDs(graph), 20);
							}


						}




						// Verify that this path can be properly reassigned in the graph.

						debugMes("\nVerifying that restructured path: " + path + " is rethreaded through the graph with fewer repeat units.", 12);
						List<Integer> updated_path = reassign_restructured_path_in_graph(graph, path);

				


						if (updated_path == null 
								|| 

								(   path.size() > 1 &&	
										(updated_path.equals(path)
												||
												(Path.countNumNodesNotUnique(updated_path) != 0)
												)

										)

								) {

							throw new RuntimeException("Remaining terminal repeat-unrolled path: " + path 
									+ " " + Path.getRepeatNodesAndCounts(path)
									+ " was not properly restructured: "
									+ updated_path + " " + Path.getRepeatNodesAndCounts(updated_path));

						}
						else {
							debugMes("\tVerification OK: path:" + path + " " + Path.getRepeatNodesAndCounts(path) + 
									" => " +  updated_path + " " + Path.getRepeatNodesAndCounts(updated_path) + "\n", 12); 
						}

					} // end if restructured_flag
					
					
				} // end for path
				
				
			} // end for pp
			
			
	
			
			if (true) {
				String filename = FILE + ".TerminalRepeatUnroll_" + unroll_loop_counter + ".dot";
				PrintStream p;
				try {
					p = new PrintStream(new FileOutputStream(filename));
					writeDotFile(graph,p,"TerminalRepeatUnroll_" + unroll_loop_counter);
					p.close();
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
			}
			
			
			// reassign repeat nodes to their new nodes in the graph
			debugMes("\n\n## Post-terminal-repeat-unroll round: " + unroll_loop_counter + ", reassigning_repeat_nodes_in_reads\n", 10);

			int num_paths_updated = reassign_repeat_nodes_in_reads(graph, combinedReadHash, 
					restructured_nodes, null, true);

			if (num_paths_updated == 0)
				throw new RuntimeException("Error, no paths were updated after this round of repeat unrolling");
			
			debugMes("\n\nNumber of paths refined: " + num_paths_updated, 10);
			
			// look for remaining repeats
			repeat_node_id_to_longest_path = find_repeat_containing_pairpaths_ignoreLastNode(combinedReadHash);
			
			
		
		}
		
		/*
		// prune the extra edges that are now not supported by the repeat-unrolled reads.
		debugMes("\n\n## Post-unroll round: " + unroll_loop_counter + ", redefine_all_graph_edges()\n", 10);
		redefine_all_graph_edges(graph, combinedReadHash);
		*/
		
		int num_paths_redefined = 1;

		HashMap<PairPath,Boolean> pp_remains_unchanged_skip_list = new HashMap<PairPath,Boolean>(); // now we use it. Graph isn't changing, so wont' expect paths to change unless they have a better placement given the unrolled edge weights.

		int refinement_round = 0;
		
		while (num_paths_redefined > 0) {

			refinement_round++;

			num_paths_redefined = reassign_repeat_nodes_in_reads(graph, 
					combinedReadHash, all_repeat_related_nodes, 
					pp_remains_unchanged_skip_list, false);

			debugMes("unroll_remaining_terminal_loops::PATH_REFINEMENT_ROUND: " + refinement_round + " NUMBER_PATHS_REDEFINED: " + num_paths_redefined, 10);

			/*
			if (num_paths_redefined > 0) {
				redefine_all_graph_edges(graph, combinedReadHash);
			}
			
			*/
		}


		return(true);
	}

	
	
	
	
	private static void redefine_all_graph_edges(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		
		// capture the edges supported by the repeat-unrolled paths
		HashMap<SimpleEdge, Integer> edge_pp_counter = new HashMap<SimpleEdge, Integer>();
		
		HashMap<SimpleEdge,String> edge_text = new HashMap<SimpleEdge,String>();
		
		
		for (HashMap<PairPath, Integer> hmap_pp : combinedReadHash.values()) {
			
			for (PairPath pp : hmap_pp.keySet()) {
				
				for (List<Integer> path : pp.get_paths()) {
					
					Integer read_support = hmap_pp.get(pp);
					
					for (int i = 0; i < path.size() - 1; i++) {
						SeqVertex from_v = getSeqVertex(graph, path.get(i));
						SeqVertex to_v = getSeqVertex(graph, path.get(i+ 1));
						SimpleEdge se = graph.findEdge(from_v,  to_v);
						if (se == null) {
							throw new RuntimeException("Error, should have edge between " + from_v + " and " + to_v + " but could not be found in graph.");
						}
						
						int edge_support = read_support;
						
						if (edge_pp_counter.containsKey(se)) {
							edge_support = edge_pp_counter.get(se) + read_support;
						}
					
							
						edge_pp_counter.put(se, edge_support);
						
						edge_text.put(se, "edge:"+from_v + "<->" + to_v + " support: " + edge_support);
						
						
					}
				}
			}
			
		}
		
		// purge all edges not found in repeat-unrolled reads
		// and reset weight of edges according to the read support of pairpaths containing them.
		ArrayList<SimpleEdge> all_edges = new ArrayList<SimpleEdge>(graph.getEdges());
		for (SimpleEdge se : all_edges) {
			if (edge_pp_counter.containsKey(se)) {
				debugMes("-Retaining edge: " + edge_text.get(se), 15);
				se.setWeight(edge_pp_counter.get(se));
			}
			else {
				debugMes("-Pruning edge: " + se, 20);
				graph.removeEdge(se);
			}
		}
		
	}
	

	private static int reassign_repeat_nodes_in_reads(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash,
			HashMap<Integer, Boolean> restructured_nodes,
			HashMap<PairPath,Boolean> pp_remains_unchanged_skip_list, // dont want it, can use null
			Boolean restrict_to_unrolled_repeat_containing_paths
			) {
	
		Set<Integer> restructured_node_ids = restructured_nodes.keySet();
		
		debugMes("Restructured nodes list: " + restructured_node_ids, 15);
		
		
		List<PairPath> orig_pps = new ArrayList<PairPath>();
		List<PairPath> updated_pps = new ArrayList<PairPath>();
		List<Integer> orig_counts = new ArrayList<Integer>();
		
		////////////////////////////////////////////////////////////////
		// iterate through all the reads and reassign where necessary.
		////////////////////////////////////////////////////////////////
		
		int pp_counter = 0;
		int total_pp = count_pp_in_combinedReadHash(combinedReadHash);
		
		for (HashMap<PairPath,Integer> pp_hmap : combinedReadHash.values()) {
			for (PairPath pp : pp_hmap.keySet()) {
				
				if (pp_remains_unchanged_skip_list != null && pp_remains_unchanged_skip_list.containsKey(pp))
					continue;
				
				
				pp_counter++;
				System.err.print("\rpp[" + pp_counter + "]  / " + total_pp + " = " + (int) ( (float)pp_counter/total_pp * 100) + " %    ");
				
				PairPath updated_pp = new PairPath(pp);
				Integer orig_count = pp_hmap.get(pp);
				
				boolean pp_updated_flag = false;
				boolean pp_not_remapped_flag = false;
				
				for (List<Integer> path : updated_pp.get_paths()) {


					if (restrict_to_unrolled_repeat_containing_paths) {
						
						if (Path.getRepeatNodesAndCounts(path).size() == 0) {
							continue; // dont try to remap it.
						}
						
					}
					
					
					if (Path.contains_any_node_id(path, restructured_node_ids)) {
						debugMes("Attempting to reassign repeat-node containing path: " + path, 15);

						List<Integer> updated_path = reassign_restructured_path_in_graph(graph, path);
						if (updated_path == null) {
							
							// shouldn't happen now.
							
							pp_not_remapped_flag = true;
							
						}
						else {
							if ( (! updated_path.equals(path)) 
									//&&
									//Path.countNumNodesNotUnique(path) > Path.countNumNodesNotUnique(updated_path) 
									&&
									score_path_by_repeats(updated_path, graph) < score_path_by_repeats(path, graph)
									
									) {

								debugMes("REASSIGNED_PATH: " + path + " " + Path.getRepeatNodesAndCounts(path)
										+ " => " + updated_path + " " + Path.getRepeatNodesAndCounts(updated_path), 15);

								// see if we just moved from a non-self terminal repeat to some other internal repeat arrangement:
								if (Path.hasTerminalNonSelfRepeat(path) && (! Path.hasTerminalNonSelfRepeat(updated_path))
										&& 
										Path.countNumNodesNotUnique(path) <= Path.countNumNodesNotUnique(updated_path) ) {

									//FIXME: better understand this extremely rare edge case
									
									// just remove the terminal repeat node, since we were unable to resolve it
									debugMes("WARNING: terminal repeat node containing path just rearranged to include alternate repeat structures that should have already been resolved earlier.", 15);
									path.remove(path.size()-1);
								}
								else {

									path.clear();
									path.addAll(updated_path);
								}
								pp_updated_flag = true;

							}
							else {
								debugMes("Path " +  path + " " + Path.getRepeatNodesAndCounts(path) 
										+ " remains unchanged or repeat count stayed the same => "
										+ updated_path + " " + Path.getRepeatNodesAndCounts(updated_path), 15);
								HashMap<Integer,Integer> rpt_nodes = Path.getRepeatNodesAndCounts(path);

								if (rpt_nodes.size() > 0) {
									debugMes("\t** path still contains repeat nodes: " + rpt_nodes, 15);
								}
							}
						}

					}
						
				}
				if (pp_updated_flag || pp_not_remapped_flag) {
					if (pp_updated_flag) {
						updated_pps.add(updated_pp);
						orig_pps.add(pp);
						orig_counts.add(orig_count);
						debugMes("PATH updated for : " + pp + " to " + updated_pp + " orig_first: " + pp.getFirstID() + ", updated_pp.first: " + updated_pp.getFirstID(), 15);
					}
					else if (pp_not_remapped_flag) {
						orig_pps.add(pp);
						updated_pps.add(null);
						orig_counts.add(-1);
						
					}
				}
				else {
					if (pp_remains_unchanged_skip_list != null) {
						pp_remains_unchanged_skip_list.put(pp, true);
					}
				}
				
			}
			
		}
		
		/////////////////////////////////////////////////////////////////////////////////////
		// reorganize any changes in the combinedReadHash based on original node identifiers
		////////////////////////////////////////////////////////////////////////////////////
		
		for (int i = 0; i < updated_pps.size(); i++) {
			PairPath updated_pp = updated_pps.get(i);
			
			// only use the orig_pp to get the first node, since data structures revolve around the actual hashmap objs.
			PairPath orig_pp = orig_pps.get(i);
			Integer orig_count = orig_counts.get(i);
			debugMes("Reorganizing combined read hash for:  orig: " + orig_pp + " to updated_pp: " + updated_pp, 15);
			
			// remove the orig pp
			Integer orig_first_node = orig_pp.getFirstID();
			
			

			combinedReadHash.get(orig_first_node).remove(orig_pp);
			if (combinedReadHash.get(orig_first_node).size() == 0) {
				combinedReadHash.remove(orig_first_node);
			}

			if (updated_pp != null) {
				// add the new pp
				Integer new_first_node = updated_pp.getFirstID();
				if (combinedReadHash.containsKey(new_first_node)) {
					combinedReadHash.get(new_first_node).put(updated_pp, orig_count);
				}
				else {
					combinedReadHash.put(new_first_node, new HashMap<PairPath,Integer>());
					combinedReadHash.get(new_first_node).put(updated_pp, orig_count);
				}
			}
			
		}
		
		return(updated_pps.size());
	}


	private static int count_pp_in_combinedReadHash(
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		
		int count = 0;
		
		for (HashMap<PairPath,Integer> pp_hmap : combinedReadHash.values()) {
			for (PairPath pp : pp_hmap.keySet()) {
				count++;
			}
		}
		
		return(count);
		
	}


	private static List<Integer> reassign_restructured_path_in_graph(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, List<Integer> path) {


		// find a complete path in the graph.
		// prefer the one with the fewest number of repeated nodes.
		
		if (path.size() == 1) {
			
			return(path); //FIXME: should try to assign it to it's best repeat node if it's a repeat, to keep it tidy.
			
		}
		
		//int repeat_cap = Path.countNumOrigNodesNotUnique(path);
		int repeat_cap = Path.countNumNodesNotUnique(path);
		
		debugMes("reassign_restructed_path_in_graph(" + path + " with cap of " + repeat_cap + " num local repeats.",15);
		
		
		
		int max_num_local_repeats = repeat_cap;

		List<Integer> chosen_thus_far = new ArrayList<Integer>();

		PATH_COUNTER = 0; // init global

		HashMap<String,List<List<Integer>>> memoize_best_path = new HashMap<String,List<List<Integer>>>();
		
		List<List<Integer>> complete_path = recursively_explore_graph_paths(graph, chosen_thus_far, path, 
				0, 0, 
				max_num_local_repeats, memoize_best_path);
		
		if (complete_path != null) {

			return(complete_path.get(0)); // take the first one of tied entries.

		}

		

		
		throw new RuntimeException("Error, couldn't remap path: " + path + " within the graph");
		
		/*
		debugMes("WARNING: couldn't remap path: " + path + " within the graph", 12);
		
		return(null); // no remapping
		 */

	}

	
	private static List<List<Integer>> recursively_explore_graph_paths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			List<Integer> chosen_thus_far, 
			List<Integer> path,
			Integer num_repeat_nodes, 
			Integer num_out_of_order, 
			int MAX_NUM_LOCAL_REPEATS, 
			HashMap<String,List<List<Integer>>> memoize_best_path) {

		// int MAX_SEARCHES_FOR_PATH_REFINEMENT = 5; //FIXME: make this a global and command-line parameter

		debugMes("recursively_explore_graph_paths(): pathLen: " + path.size() + ", chosen thus far: " + chosen_thus_far, 20);


		String curr_node_pos_token = null;
		if (chosen_thus_far.size() > 0) {
			curr_node_pos_token = "" + chosen_thus_far.size() + "_" + chosen_thus_far.get(chosen_thus_far.size()-1);


			if (memoize_best_path.containsKey(curr_node_pos_token)) {
				return(Path.clone(memoize_best_path.get(curr_node_pos_token)));
			}


		}

		boolean local_debug = false;

		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 25)
			local_debug = true;

		if (num_repeat_nodes > MAX_NUM_LOCAL_REPEATS) {
			debugMes("\t** terminating extension, max num local repeats encountered: " + num_repeat_nodes, 20);
			return (null);
		}

		List<List<Integer>>min_repeat_reconstructed_path_list = new ArrayList<List<Integer>>();
		Float min_repeat_reconstructed_path_repeat_score = null;  


		List<List<Integer>>all_possible_path_reconstructions_seen = new ArrayList<List<Integer>>(); // for debugging purposes.
		
		

		if (local_debug) {
			System.err.println("RECURSIVELY_EXPLORE_GRAPH_PATHS: chosen_thus_far: " + chosen_thus_far + ", path: " + path);
		}

		if (chosen_thus_far.size() == path.size()) {
			debugMes("\trecursion base case, found path: " + path, 20);

			// return empty list of paths to signal base case.
			return(Path.create_empty_path_list());
			
			//return(Path.create_list_of_paths_from_single_node_id(chosen_thus_far.get(chosen_thus_far.size()-1))); // done, return last node.
		}


		Integer target_path_node_id = path.get(chosen_thus_far.size());
		SeqVertex target_path_node = getSeqVertex(graph, target_path_node_id);
		Integer current_orig_node_id = target_path_node.getOrigButterflyID();


		// get the last node chosen thus far
		SeqVertex last_node = null;
		if (chosen_thus_far.size() > 0) {
			last_node = getSeqVertex(graph, chosen_thus_far.get(chosen_thus_far.size()-1));
			debugMes("EXTENDING FROM LAST_NODE: " + last_node + ", searching for an origID: " + current_orig_node_id, 20);
		}


		//List<SeqVertex> candidate_vertices = SeqVertex.getAllNodesHavingOriginalID(current_orig_node_id);
		
		// only pursue multiple candidates if it's an unrolled repeat vertex in the target path
		// if it's a repeat vertex that was already unrolled, add the unrolled node here.
		List<SeqVertex> candidate_vertices = SeqVertex.getAllNodesHavingOriginalID(target_path_node_id);
		if (candidate_vertices.isEmpty()) {
			candidate_vertices.add(target_path_node);
		}
		
		
		if (last_node == null) {

			/*
			if (candidate_vertices.size() > MAX_SEARCHES_FOR_PATH_REFINEMENT) {
				debugMes("Not seeding on repetitive node, skipping this path: " + path, 12); //FIXME: instead, redo seeding on non-repetitive node of this path.
				return(null);
			}
			 */


			debugMes("Initial candidate vertices based on orig_id: " + current_orig_node_id + " are " +  candidate_vertices, 20); 


		}

		else {
			// restrict the list to those vertices that are direct successors to last node
			// ensure there's an edge between previously chosen node and this one.

			candidate_vertices = last_node.getListOfSuccessors(graph, candidate_vertices);

			debugMes("\tFiltered candidate vertices for extension from: " + last_node.getID() + " are " +  candidate_vertices, 20); 



			// now sort them by unrolled edge weights.
			//if (candidate_vertices.size() > MAX_SEARCHES_FOR_PATH_REFINEMENT) {

			HashMap<Integer,Double> unrolled_weights = new HashMap<Integer,Double>();
			for (SeqVertex v : candidate_vertices) {
				SimpleEdge se = graph.findEdge(last_node, v);
				double unrolled_weight = se.get_repeat_unroll_weight();
				unrolled_weights.put(v.getID(), unrolled_weight);


			}
			final HashMap<Integer,Double> unrolled_weights_map = new HashMap<Integer,Double>(unrolled_weights);

			Collections.sort(candidate_vertices, new Comparator<SeqVertex>() {

				public int compare (SeqVertex v_a, SeqVertex v_b) {

					// want most highly supported pairpaths to sort descendingly

					Integer v_a_id = v_a.getID();
					Integer v_b_id = v_b.getID();

					if (unrolled_weights_map.get(v_a_id) < unrolled_weights_map.get(v_b_id)) {
						return(-1);
					}
					else if (unrolled_weights_map.get(v_a_id) > unrolled_weights_map.get(v_b_id)) {
						return(1);
					}
					else {
						return(0);
					}
				}

			});


			//candidate_vertices = candidate_vertices.subList(0, MAX_SEARCHES_FOR_PATH_REFINEMENT);

			//debugMes("-restricting recursive search from " + last_node.getID() + " to " + candidate_vertices, 15);



			//}


		}


		if (local_debug) {
			System.err.println("RECURSIVELY_EXPLORE_GRAPH_PATHS: candidate_vertices with orig_node_id: " 
					+ current_orig_node_id + " are: " + candidate_vertices);
		}


		int search_refinement_count = 0;



		for (SeqVertex v : candidate_vertices) {


			int local_num_repeats = num_repeat_nodes;
			int local_num_out_of_order = num_out_of_order;

			// all repeats should have been expanded by now.
			if (chosen_thus_far.contains(v.getID())) {
				if (local_debug) {
					System.err.println("already chose id: " + v.getID());
				}

				local_num_repeats++;

			}
			if (last_node != null && v.getNodeDepth() < last_node.getNodeDepth()) {
				local_num_out_of_order++;
			}


			////////////////////////////////////////////////
			// Take this node and explore other extensions.
			////////////////////////////////////////////////

			search_refinement_count++;

			/*
			if (search_refinement_count > MAX_SEARCHES_FOR_PATH_REFINEMENT)
				break;
			 */

			chosen_thus_far.add(v.getID());
			List<List<Integer>> tied_reconstructed_paths = recursively_explore_graph_paths(graph, chosen_thus_far, path, 
					local_num_repeats, local_num_out_of_order, MAX_NUM_LOCAL_REPEATS,
					memoize_best_path);

			if (tied_reconstructed_paths != null) {


				// include current node in path before scoring

				tied_reconstructed_paths = Path.prepend_node_id_to_paths(chosen_thus_far.get(chosen_thus_far.size()-1), tied_reconstructed_paths); // add current node to the lowest repeat extension.

				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
					all_possible_path_reconstructions_seen.addAll(tied_reconstructed_paths); // for debugging
				}
				debugMes("\nAll Paths and scores:", 15);
				
				for (List<Integer> reconstructed_path : tied_reconstructed_paths) {

					float repeated_node_score = score_path_by_repeats(reconstructed_path, graph);

					debugMes("score:" + repeated_node_score + " " + reconstructed_path + " " + Path.getRepeatNodesAndCounts(reconstructed_path), 15);
					
					
					if (min_repeat_reconstructed_path_repeat_score == null 
							|| Math.abs(min_repeat_reconstructed_path_repeat_score - repeated_node_score) < 0.00001)  // consider a tie
					{
						if (min_repeat_reconstructed_path_repeat_score == null || repeated_node_score < min_repeat_reconstructed_path_repeat_score ) {
							min_repeat_reconstructed_path_repeat_score = repeated_node_score;
						}

						min_repeat_reconstructed_path_list.add(reconstructed_path);

					}
					else if (repeated_node_score < min_repeat_reconstructed_path_repeat_score) {

						min_repeat_reconstructed_path_list.clear(); // reset since have lower score (better)
						min_repeat_reconstructed_path_list.add(reconstructed_path);
						min_repeat_reconstructed_path_repeat_score = repeated_node_score;
					}


				}

				
			}

			chosen_thus_far.remove(chosen_thus_far.size()-1); // remove the last element added in prep for the next one.
			
		}
		
		
		if (! min_repeat_reconstructed_path_list.isEmpty()) {

			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
				debugMes("\nALL CANDIDATE PATHS SEEN AT " + curr_node_pos_token + ":", 15);
				if (all_possible_path_reconstructions_seen.size() > 1) {
					debugMes("MULTIPLE CANDIDATE PATHS SEEN AT NODE", 15);
				}
				for (List<Integer> candidate_path : all_possible_path_reconstructions_seen) {
					float candidate_path_score = score_path_by_repeats(candidate_path, graph);
					debugMes("score: " + candidate_path_score + " "  + candidate_path + " " + Path.getRepeatNodesAndCounts(candidate_path), 15);
				}
				
			}
			
			
			debugMes("\nMinRepeat tied paths of length: " + min_repeat_reconstructed_path_list.get(0).size() + " with score: " + min_repeat_reconstructed_path_repeat_score + ":", 15);
			for (List<Integer> reconstructed_path : min_repeat_reconstructed_path_list) {
				debugMes(reconstructed_path + " " + Path.getRepeatNodesAndCounts(reconstructed_path), 15);
			}

			memoize_best_path.put(curr_node_pos_token, Path.clone(min_repeat_reconstructed_path_list));

			return(min_repeat_reconstructed_path_list);

		}

		return(null); // no min repeat paths to report.

	}
	

	private static float score_path_by_repeats(
			List<Integer> reconstructed_path,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		float OUT_OF_ORDER_PENALTY_FACTOR = 1000;
		float UNROLLED_EDGE_USE_SCORE_FACTOR = 100;
		
		// found a good path.
		
		float repeated_node_score = Path.countNumNodesNotUnique(reconstructed_path);	 
		
		
		// further penalize by the number of out-of-order nodes
		repeated_node_score += Path.countNumOutOfOrder(graph, reconstructed_path) / OUT_OF_ORDER_PENALTY_FACTOR;
		
		// take into account use of repeat-unrolled edges
		double unrolled_edge_weight_sum = Path.getUnrolledEdgeWeightSum(graph, reconstructed_path);
		repeated_node_score -= unrolled_edge_weight_sum / UNROLLED_EDGE_USE_SCORE_FACTOR;
		
		return(repeated_node_score);
		
		
	}


	private static HashMap<Integer, PairPath> find_repeat_containing_pairpaths (
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		HashMap<Integer, PairPath> repeat_node_id_to_longest_path = new HashMap<Integer, PairPath>();
		HashMap<Integer, Integer> repeat_node_id_to_max_repeat_count = new HashMap<Integer,Integer>();
		
		for (HashMap<PairPath, Integer> pp_n_counts : combinedReadHash.values()) {
			
			for (PairPath pp : pp_n_counts.keySet()) {
				
				// get all pp containing a repeat node.
				HashMap<Integer, Integer> repeat_node_ids_n_counts = pp.getRepeatNodesAndCounts();
				
				
				// assign each repeat node to the path that contains it as a repeat the greatest number of occurrences.
				if (repeat_node_ids_n_counts.size() > 0) {
					
					debugMes("repeat_node_ids_n_counts: " + repeat_node_ids_n_counts + " , pp: " + pp + ", counts: " + pp_n_counts.get(pp), 14);
					
					
					for (Integer node_id : repeat_node_ids_n_counts.keySet()) {
						Integer count = repeat_node_ids_n_counts.get(node_id);
						if ( (! repeat_node_id_to_max_repeat_count.containsKey(node_id)) 
								|| 
								repeat_node_id_to_max_repeat_count.get(node_id) < count) {
							
							repeat_node_id_to_longest_path.put(node_id, pp);
							repeat_node_id_to_max_repeat_count.put(node_id, count);
							
						}
					}
					
				}
				
			}
		}
		
		
		return(repeat_node_id_to_longest_path);
		
		
	}

	
	private static HashMap<Integer, PairPath> find_repeat_containing_pairpaths_ignoreLastNode (
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		HashMap<Integer, PairPath> repeat_node_id_to_longest_path = new HashMap<Integer, PairPath>();
		HashMap<Integer, Integer> repeat_node_id_to_max_repeat_count = new HashMap<Integer,Integer>();
		
		for (HashMap<PairPath, Integer> pp_n_counts : combinedReadHash.values()) {
			
			for (PairPath pp : pp_n_counts.keySet()) {
				
				// get all pp containing a repeat node.
				HashMap<Integer, Integer> repeat_node_ids_n_counts = pp.getRepeatNodesAndCounts_ignoreLastNode();
				
				
				// assign each repeat node to the path that contains it as a repeat the greatest number of occurrences.
				if (repeat_node_ids_n_counts.size() > 0) {
					
					debugMes("repeat_node_ids_n_counts: " + repeat_node_ids_n_counts + " , pp: " + pp + ", counts: " + pp_n_counts.get(pp), 14);
					
					
					for (Integer node_id : repeat_node_ids_n_counts.keySet()) {
						Integer count = repeat_node_ids_n_counts.get(node_id);
						if ( (! repeat_node_id_to_max_repeat_count.containsKey(node_id)) 
								|| 
								repeat_node_id_to_max_repeat_count.get(node_id) < count) {
							
							repeat_node_id_to_longest_path.put(node_id, pp);
							repeat_node_id_to_max_repeat_count.put(node_id, count);
							
						}
					}
					
				}
				
			}
		}
		
		
		return(repeat_node_id_to_longest_path);
		
		
	}
	

	private static void infer_best_triplets_across_unresolved_Xstructure(
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, Boolean> xStructuresResolvedByTriplets,
			HashMap<Integer, List<List<Integer>>> tripletMapper) {
	
			for (Integer xstructure_node : xStructuresResolvedByTriplets.keySet()) {
				if (! xStructuresResolvedByTriplets.get(xstructure_node)) {
					
					debugMes("Examining unresolved X structure at: " + xstructure_node, 10);
					
					SeqVertex v = getSeqVertex(graph, xstructure_node);
					
					Iterator<SeqVertex> predecessors = graph.getPredecessors(v).iterator();
					Iterator<SeqVertex> successors = graph.getSuccessors(v).iterator();
					
					
					
					while(predecessors.hasNext()) {
						SeqVertex p = predecessors.next();
					
						if (successors.hasNext()) {
							SeqVertex s = successors.next();
							
							List<Integer> triplet = new ArrayList();
							triplet.add(p.getID());
							triplet.add(xstructure_node);
							triplet.add(s.getID());
							
							if (! tripletMapper.containsKey(xstructure_node))
								tripletMapper.put(xstructure_node, new ArrayList<List<Integer>>());
							
							tripletMapper.get(xstructure_node).add(triplet);
							
							debugMes("INFERRING triplet for UNRESOLVED X STRUCTURE (" + xstructure_node + ") -> " + triplet, 10);
					
						
						}
					}
					
				}
			}
		
		
	}

	/*

	private static void printGraph (
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		for (SeqVertex v : graph.getVertices()) {
			System.out.println("Vertex: " + v.getID() + ", seq: " + v.getName());
		}
		
	}
    */



	private static void reduce_to_max_paths_per_node(
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash,
			int max_num_paths_per_start_node) {
		

		for (Integer start_node : componentReadHash.keySet() ) {

			final HashMap<PairPath,Integer> pp_to_counts = componentReadHash.get(start_node);
			
			List<PairPath> pair_paths_list = new ArrayList<PairPath>(pp_to_counts.keySet());
			
			if (pair_paths_list.size() > max_num_paths_per_start_node) {
				
			
				Collections.sort(pair_paths_list, new Comparator<PairPath>() {

					public int compare (PairPath pp_A, PairPath pp_B) {

						// want most highly supported pairpaths to sort descendingly

						int read_support_A = pp_to_counts.get(pp_A);
						int read_support_B = pp_to_counts.get(pp_B);

						if (read_support_A > read_support_B)
							return -1;
						else if (read_support_A < read_support_B)
							return 1;
						else
							return 0;

					}

				});

				List<PairPath> to_remove = pair_paths_list.subList(max_num_paths_per_start_node, pair_paths_list.size());
				
				for (PairPath pp : to_remove)
					componentReadHash.get(start_node).remove(pp);
			
				
			}
			
			
			

		}
		
		
	}


	private static void report_pairpath_counts(
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash) {
		
		
		System.out.println("*****  PairPath Counts *****");
		
		int count_of_total_pps = 0;
		
		for (Integer start_node : componentReadHash.keySet() ) {
			
			int count_pp_at_node = 0;
			
			String indiv_read_support_text = "";
			
			final HashMap<PairPath,Integer> pp_to_counts = componentReadHash.get(start_node);
			
			debugMes("componentReadHash, start node: " + start_node + " has size: " + pp_to_counts.size(), 12);
			
			List<PairPath> pair_paths_list = new ArrayList<PairPath>(pp_to_counts.keySet());
			for (PairPath pp : pair_paths_list) {
				
				//debugMes("CHECKING-A: " + pp, 12);
				
				Integer read_support = componentReadHash.get(start_node).get(pp);
				
				//debugMes("CHECKING-B: " + pp + " has read support: " + read_support, 12);
				
				//debugMes("pp: " + pp + " has read support: " + read_support, 12);
				
				if (read_support == null) {
					componentReadHash.get(start_node).put(pp, 1); //FIXME: shouldn't have null entries here.
					debugMes("\tERROR: no support for pp: " + pp, 12);
				}
			}
			
			Collections.sort(pair_paths_list, new Comparator<PairPath>() {

				public int compare (PairPath pp_A, PairPath pp_B) {

					// want most highly supported pairpaths to sort descendingly

					int read_support_A = pp_to_counts.get(pp_A);
					int read_support_B = pp_to_counts.get(pp_B);

					if (read_support_A > read_support_B)
						return -1;
					else if (read_support_A < read_support_B)
						return 1;
					else
						return 0;

				}

			});

			
			
			
			for (PairPath pp : pair_paths_list) {
				
				int read_support = componentReadHash.get(start_node).get(pp);
				
				indiv_read_support_text += "\t" + pp + " has read support: " + read_support + "\n";
				
				count_pp_at_node++;
			}
			
			System.out.println("Node: " + start_node + " has " + count_pp_at_node + " pairpaths stored:\n" + indiv_read_support_text);
			
			count_of_total_pps += count_pp_at_node;
			
		}
		
		System.out.println("## Total number of pairpaths: " + count_of_total_pps);
		
	}


	private static int handleRemainingCyclicReads(
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash, DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
	
		
		int count_of_fractured_reads = 0;
		
		
		HashMap<PairPath,Integer> readParts = new HashMap<PairPath,Integer>();
		
		for (Integer start_node : componentReadHash.keySet() ) {
			
			List<PairPath> to_purge = new ArrayList<PairPath>();
			for (PairPath pp : componentReadHash.get(start_node).keySet()) {
		
				// check for fractured path
				boolean fractured = false;
	
				for (List<Integer> path : pp.get_paths()) {
					Iterator<Integer> i = path.iterator();
				
					int prev_depth = -1;
					HashMap<Integer,Boolean> node_visitor = new HashMap<Integer,Boolean>();
					while (i.hasNext()) {
						Integer node_id = i.next();
						SeqVertex v = getSeqVertex(graph, node_id);
						Integer depth = v.getNodeDepth();
						if (prev_depth > depth || node_visitor.containsKey(node_id)) {
							fractured = true;
							break;
						}
						
						prev_depth = depth;
						node_visitor.put(node_id, true);
						
					}
					if (fractured)
						break;
				}
				if (fractured) {
					to_purge.add(pp);
				}
				else {
					// double check there's no cycle here:
					if (pp.max_count_occurrence_individual_node_in_path(pp) > 1) {
						throw new RuntimeException("Error, path:" + pp + " involves an undetected cycle");
					}
				}
			}
			
			if (! to_purge.isEmpty()) {
				
				for (PairPath pp : to_purge) {
					count_of_fractured_reads++;
					debugMes("DAG_CONFLICTING_READ_FRAGMENTED: " + pp, 10);
					componentReadHash.get(start_node).remove(pp);
					List<List<Integer>> parts = fragment_DAG_conflicting_pairpath(pp, graph);
					for (List<Integer> read_part : parts) {
						
						PairPath part_pairpath = new PairPath(read_part);
						if (readParts.containsKey(part_pairpath)) {
							readParts.put(part_pairpath, readParts.get(part_pairpath) + 1);
						}
						else {
							readParts.put(part_pairpath, new Integer(1));
						}
					}
					
				}
			}
		
		}
		
		// add the fragments back in.
		for (PairPath pp : readParts.keySet()) {
			Integer first_id = pp.getFirstID();
			if (! componentReadHash.containsKey(first_id)) {
				componentReadHash.put(first_id, new HashMap<PairPath,Integer>());
			}
			if (! componentReadHash.get(first_id).containsKey(pp)) {
				componentReadHash.get(first_id).put(pp, 1);
			}
			else {
				componentReadHash.get(first_id).put(pp, componentReadHash.get(first_id).get(pp)+1);
			}
		}
	
		
		return count_of_fractured_reads;
	}


	private static List<List<Integer>> fragment_DAG_conflicting_pairpath(
			PairPath pp, DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		List<List<Integer>> read_parts = new ArrayList<List<Integer>>();
		
		List<Integer> part = new ArrayList<Integer>();
		
		List<List<Integer>> node_depths_tracker = new ArrayList<List<Integer>>();
		for (List<Integer> path : pp.get_paths()) {
		
			if (path.isEmpty())
				continue;
			
			ArrayList<Integer> node_depths_list = new ArrayList<Integer>();
			
			Iterator<Integer> it = path.iterator();
			int prev_depth = -1;
			
			HashMap<Integer,Boolean> node_visitor = new HashMap<Integer,Boolean>();
			while (it.hasNext()) {
				Integer node_id = it.next();
				
				SeqVertex v = getSeqVertex(graph, node_id);
				if (v._node_depth < prev_depth || node_visitor.containsKey(node_id)) {
					// problem...
					// fracture here.
					if (! part.isEmpty()) {
						read_parts.add(part);
						part = new ArrayList<Integer>();
						
					}
				}
				node_depths_list.add(v._node_depth);
				prev_depth = v._node_depth;
				node_visitor.put(node_id, true);
				part.add(node_id);
			}
			if (! part.isEmpty()) {
				read_parts.add(part);
				part = new ArrayList<Integer>();
			}
			
			node_depths_tracker.add(node_depths_list);
		}
		
		debugMes("FRACTURED pairpath: " + pp + " with node_depths: " + node_depths_tracker  + " into " + read_parts.size() + ": " + read_parts, 10);
		
		return read_parts;
	}
	

	private static HashMap<List<Integer>, Pair<Integer>> pasafly(
			final DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis, 
			HashMap<Integer, List<List<Integer>>> tripletMapper, HashMap<Integer, List<List<Integer>>> extendedTripletMapper) {
	
		
		
		debugMes("Beginning PasaFly",10);
	
		PasaVertex.max_top_paths_to_store = TransAssembly_allProbPaths.MAX_NUM_PATHS_PER_PASA_NODE;
		
		// populate pairPathToReadSupport: PairPath => readSupport
		// and pairPaths hashset: the list of all PairPaths
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(componentReadHash, pairPaths, pairPathToReadSupport);
			
		ArrayList<PairPath> pairPathsSortedList = new ArrayList<PairPath>(pairPaths);
		
		
		
		
		Comparator<PairPath> pairPathOrderComparer = new Comparator<PairPath>() { // sort by first node depth in graph
			public int compare(PairPath a, PairPath b) {
				
				
				if (a.equals(b)) {
					return(0);
				}
				
				// check first node
				
				//----------------------
				// use node depth in graph
				//----------------------
				
				// check first node

			Integer a_index = a.getFirstID();
			Integer b_index = b.getFirstID();

			int f1 = getSeqVertex(graph, a_index)._node_depth;   // why using FinishingTime instead of DiscoveryTime?
			int f2 = getSeqVertex(graph, b_index)._node_depth;
			if( f1 < f2 )
				return -1;
			else if( f1 > f2 )
				return 1;

			// same node depth.
			if (a_index < b_index)
				return -1;
			else if (a_index > b_index)
				return 1;
			
			// if here, 
			// same first node ID
			
		
			
			// check last node
			Integer a_last_index = a.getLastID();
			Integer b_last_index = b.getLastID();

			int l1 = getSeqVertex(graph,a_last_index)._node_depth;
			int l2 = getSeqVertex(graph,b_last_index)._node_depth;
			if (l1 < l2) {
				return(-1);
			}
			else if (l1 > l2) {
				return(1);
			}
			
			// same last node depth too.
			// compare their node identifiers
			if (a_last_index < b_last_index)
				return(-1);
			else if (a_last_index > b_last_index)
				return(1);
			
			
			// default
			// not the same paths, but same start node and last node DFS, so just order based on hashcode
			return ( (a.hashCode() < b.hashCode()) ? 1 : -1);
				
		}
	
	};
	
				
				
		
		
		Collections.sort(pairPathsSortedList, pairPathOrderComparer);
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			debugMes("SORTED PAIRPATHS IN ORDER:", 15);
			for (PairPath p : pairPathsSortedList) {
				debugMes("\t" + p, 15);
			}
			
			
		}
		
		
		ArrayList<PasaVertex> pasaVerticesSortedList = new ArrayList<PasaVertex>();
		for (PairPath pp : pairPathsSortedList) { // already sorted
			int count = pairPathToReadSupport.get(pp);
			pasaVerticesSortedList.add(new PasaVertex(pp, count));
		}
		
		PasaVertex [] pasaVerticesSortedArr = pasaVerticesSortedList.toArray(new PasaVertex[pasaVerticesSortedList.size()]);
		PairPath[] pairPathsSortedArr = pairPathsSortedList.toArray(new PairPath[pairPathsSortedList.size()]);
		
		
		//----------------------
		// EXAMINE CONTAINMENTS
		//----------------------
		
		// init
		ArrayList<PairPath> pairPathsContainmentsRemoved = new ArrayList<PairPath>(pairPathsSortedList);
		ArrayList<PasaVertex> pasaVerticesContainmentsRemoved = new ArrayList<PasaVertex>(pasaVerticesSortedList);

		
		debugMes("Assigning pairpath containments.", 10);
		List<Integer> containments = assignPasaPairPathContainments(graph, dijkstraDis, pasaVerticesSortedArr); // vertices updated to include containment info.

		debugMes("REMOVING CONTAINMENTS: " + containments, 10);
		for(int i = 0; i < containments.size(); i++)
		{

			pasaVerticesContainmentsRemoved.remove(pasaVerticesSortedArr[containments.get(i)]);
			pairPathsContainmentsRemoved.remove(pairPathsSortedArr[containments.get(i)]);
		}
		
		//------------------------------------------------
		// EXAMINE UNCERTAINTIES THAT BREAK TRANSITIVITY
		//------------------------------------------------

		HashSet<Integer> vertices = extract_vertex_list_from_PairPaths(pairPathsContainmentsRemoved);
		
		PairPath [] pairPathsContainmentsRemovedArr = pairPathsContainmentsRemoved.toArray(new PairPath[pairPathsContainmentsRemoved.size()]);
		PasaVertex [] pasaVerticesContainmentsRemovedArr = pasaVerticesContainmentsRemoved.toArray(new PasaVertex[pasaVerticesContainmentsRemoved.size()]);
		
		boolean[][] dag = getPairPathConsistencyDAG(graph, dijkstraDis, pairPathsContainmentsRemovedArr);

		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 10) {
			debugMes("PASA Consistency DAG--------------------------\n", 10);

			System.out.println(boolean_matrix_toString(dag));
		}

		

		ArrayList<PairPath> pairPathsUncertainRemoved = new ArrayList<PairPath>(pairPathsContainmentsRemoved);
		ArrayList<PasaVertex> pasaVerticesUncertainRemoved = new ArrayList<PasaVertex>(pasaVerticesContainmentsRemoved);

	
		debugMes("Identifying uncertain entries that break transitivities.", 10);
		// identify and remove uncertain entries (those that break transitive compatibility relationships)
		ArrayList<Integer> uncertain = getUncertainRequireOverlap(dag, pairPathsContainmentsRemovedArr, graph, dijkstraDis);

		debugMes("Uncertain indices include: " + uncertain, 10);

		
		debugMes("REMOVING UNCERTAINTIES: " + uncertain, 10);

		for(int i = 0; i < uncertain.size(); i++)
		{

			pasaVerticesUncertainRemoved.remove(pasaVerticesContainmentsRemovedArr[uncertain.get(i)]);
			pairPathsUncertainRemoved.remove(pairPathsContainmentsRemovedArr[uncertain.get(i)]);
		}

		HashSet<Integer> vertices_after_removed_uncertainties = extract_vertex_list_from_PairPaths(pairPathsUncertainRemoved);
		
		if (vertices_after_removed_uncertainties.size() < vertices.size()) {
			int missing_node_count = vertices.size() - vertices_after_removed_uncertainties.size();
			debugMes("WARNING, MISSING: " + missing_node_count + " of " + vertices.size() + " nodes after removing uncertainties", 10);
			for (Integer v : vertices) {
				if (! vertices_after_removed_uncertainties.contains(v)) {
					debugMes("WARNING, MISSING NODE: After removing uncertainties, missing node from graph: " + v, 10);
				}
			}
		}
			
		
		PasaVertex[] pasaVerticesUncertainRemovedArr = pasaVerticesUncertainRemoved.toArray(new PasaVertex[pasaVerticesUncertainRemoved.size()]);
		PairPath[] pairPathsUncertainRemovedArr = pairPathsUncertainRemoved.toArray(new PairPath[pairPathsUncertainRemoved.size()]);

		//print pair paths
		debugMes("PAIR PATHS remaining after uncertainties removed ---------------------------------------",10);
		for(int i = 0; i < pairPathsUncertainRemovedArr.length; i++)
		{
			debugMes("PairPathAfterUncertainRemoved "+ i + " " + pairPathsUncertainRemovedArr[i], 10);
		}

		// regenerate the dag now that the uncertain entries are removed.
		dag = getPairPathConsistencyDAG(graph, dijkstraDis, pairPathsUncertainRemovedArr); // already identified containments

		//print dag
		debugMes("DAG after uncertainties removed ---------------------------\n" + boolean_matrix_toString(dag),10);


		//2.2 check transitivity
		if(!checkTransitivityRequireOverlap(dag, pairPathsUncertainRemovedArr, graph, dijkstraDis))
		{
			throw(new RuntimeException("Graph is NOT transitive!"));

		}
		else {
			debugMes("Transitivity of compatibility graph validates.", 10);
		}

		
		
		//2.2 check transitivity
		if(!checkTransitivityRequireOverlap(dag, pairPathsUncertainRemovedArr, graph, dijkstraDis))
		{
			throw(new RuntimeException("Graph is NOT transitive!"));
			
		}
		else {
			debugMes("Transitivity of compatibility graph validates.", 10);
		}
		
		// track the final vertex identifiers
		HashMap<PairPath,Integer> finalVertexPositions = new HashMap<PairPath,Integer>();
		for (int i = 0; i < pairPathsUncertainRemovedArr.length; i++) {
			finalVertexPositions.put(pairPathsUncertainRemovedArr[i], i);
		}
		
		
		debugMes("build_PASA_trellis_left_to_right()", 10);
		
		build_PASA_trellis_left_to_right(pasaVerticesUncertainRemovedArr, dag, graph, componentReadHash, dijkstraDis, 
					pairPathToReadSupport, tripletMapper, extendedTripletMapper);
		
		// get highest scoring path:
		debugMes("Identifying highest scoring PASA path.", 10);
		ScoredPath best = null;
		for (int i = 0; i < pasaVerticesUncertainRemovedArr.length; i++) {
			
			ScoredPath sp = pasaVerticesUncertainRemovedArr[i].get_highest_scoring_fromPath();
			if (best == null || sp.score > best.score) {
				best = sp;
			}
		}
		
		debugMes("Best score: " + best.score + ", containing entries: " + best.paths, 10);
		
		List<Integer> best_path_vertex_list = Path.collapse_compatible_pair_paths(best.paths);
		
		HashMap<List<Integer>, Pair<Integer>> final_transcripts = new HashMap<List<Integer>, Pair<Integer>>();
		
		final_transcripts.put(best_path_vertex_list, new Pair<Integer>(1,1));
		
		// remove those pairpaths included in the best path
		List<PairPath> toRemove = new ArrayList<PairPath>();
		for (PairPath pp : finalVertexPositions.keySet()) {
			
			if (pp.isCompatibleAndContainedBySinglePath(best_path_vertex_list))
				toRemove.add(pp);
		}
		for (PairPath pp : toRemove)
			finalVertexPositions.remove(pp);
		toRemove.clear();
		
		
		
		// Now, extract the top combined path that contains each missing transcript
		// Prioritize according to paired path support, and break ties according to representing the most additional missing entries.
		debugMes("build_PASA_trellis_right_to_left()", 10);
		
		build_PASA_trellis_right_to_left(pasaVerticesUncertainRemovedArr, dag, graph, componentReadHash, dijkstraDis, pairPathToReadSupport, tripletMapper);
		
		List<PairPath> unrepresented_pairpaths = new ArrayList<PairPath>(finalVertexPositions.keySet());
		
		final HashMap<PairPath,Integer>pairPathToReadSupportFixed  = new HashMap<PairPath,Integer>(pairPathToReadSupport);
		
		Collections.sort(unrepresented_pairpaths, new Comparator<PairPath>() {
			public int compare(PairPath a, PairPath b) {
				
				int count_a = pairPathToReadSupportFixed.get(a);
				int count_b = pairPathToReadSupportFixed.get(b);
			
				if (count_a == count_b) {
					return(0);
				}
				else if (count_a > count_b) {
					return(-1);
				}
				else {
					return(1);
				}
				
			}
			
		});
		
		
		Iterator<PairPath> it = unrepresented_pairpaths.iterator();
		
		while (it.hasNext() && ! finalVertexPositions.isEmpty()) {
			PairPath pp = it.next();
			
			if (! finalVertexPositions.containsKey(pp)) {
				// recovered in a previous round
				continue;
			}
		
			debugMes("Nucleating next pasa path on PP: " + pp + ", having read support: " + pairPathToReadSupportFixed.get(pp), 10);
			
			// get the highest scoring chain that contains pp
			int index = finalVertexPositions.get(pp);
			
			List<ScoredPath> sp_from = pasaVerticesUncertainRemovedArr[index].get_all_highest_scoring_fromPath();
			List<ScoredPath> sp_to = pasaVerticesUncertainRemovedArr[index].get_all_highest_scoring_toPath();
			
			debugMes("Best combined partial paths containing pairpath: " + pp + " include: (From): " 
					+ sp_from + ", (To): " + sp_to, 10);
			
			List<PairPath> combined_pp_list = new ArrayList<PairPath>();
			if (sp_from.size() > 1 || sp_to.size() > 1) {
				// find the combination that covers the most currently unrepresented pairpaths
				combined_pp_list = find_paired_paths_with_greatest_map_support(sp_from, sp_to, finalVertexPositions);
			}
			else {
				// single path each.
				combined_pp_list.addAll(sp_from.get(0).paths);
				combined_pp_list.addAll(sp_to.get(0).paths);
			}
			
			List<Integer> combined_path_vertex_list = Path.collapse_compatible_pair_paths(combined_pp_list);
			final_transcripts.put(combined_path_vertex_list, new Pair<Integer>(1,1));
			// remove those pairpaths included in the best path
		
			for (PairPath p : finalVertexPositions.keySet()) {
				
				if (p.isCompatibleAndContainedBySinglePath(combined_path_vertex_list))
					toRemove.add(p);
			}
			for (PairPath p : toRemove)
				finalVertexPositions.remove(p);
			
			toRemove.clear();
			
			
		}
		
		return(final_transcripts);
		
		
	}

	private static HashMap<List<Integer>, Pair<Integer>> pasaflyunique (
			final DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis, HashMap<Integer, List<List<Integer>>> tripletMapper, HashMap<Integer, List<List<Integer>>> extendedTripletMapper) {
	
		debugMes("Beginning PasaFlyUnique",10);
	
		PasaVertex.max_top_paths_to_store = 1; //TransAssembly_allProbPaths.MAX_NUM_PATHS_PER_PASA_NODE;
		
		// populate pairPathToReadSupport: PairPath => readSupport
		// and pairPaths hashset: the list of all PairPaths
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
		populate_pairpaths_and_readsupport(componentReadHash, pairPaths, pairPathToReadSupport);
			
		ArrayList<PairPath> pairPathsSortedList = new ArrayList<PairPath>(pairPaths);
		
		
		
		
		Comparator<PairPath> pairPathOrderComparer = new Comparator<PairPath>() { // sort by first node depth in graph
			public int compare(PairPath a, PairPath b) {
				
				
				if (a.equals(b)) {
					return(0);
				}
				
				// check first node
				
				//----------------------
				// use node depth in graph
				//----------------------
				
				// check first node

			Integer a_index = a.getFirstID();
			Integer b_index = b.getFirstID();

			int f1 = getSeqVertex(graph, a_index)._node_depth;   // why using FinishingTime instead of DiscoveryTime?
			int f2 = getSeqVertex(graph, b_index)._node_depth;
			if( f1 < f2 )
				return -1;
			else if( f1 > f2 )
				return 1;

			// same node depth.
			if (a_index < b_index)
				return -1;
			else if (a_index > b_index)
				return 1;
			
			// same first node ID
			
		
			
			// check last node
			Integer a_last_index = a.getLastID();
			Integer b_last_index = b.getLastID();

			int l1 = getSeqVertex(graph,a_last_index)._node_depth;
			int l2 = getSeqVertex(graph,b_last_index)._node_depth;
			if (l1 < l2) {
				return(-1);
			}
			else if (l1 > l2) {
				return(1);
			}
			
			// same last node depth too.
			// compare their node identifiers
			if (a_last_index < b_last_index)
				return(-1);
			else if (a_last_index > b_last_index)
				return(1);
			
			
			// default
			// not the same paths, but same start node and last node DFS, so just order based on hashcode
			return ( (a.hashCode() < b.hashCode()) ? 1 : -1);
				
		}
	
	};
	
				
				
		
		
		Collections.sort(pairPathsSortedList, pairPathOrderComparer);
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			debugMes("SORTED PAIRPATHS IN ORDER:", 15);
			for (PairPath p : pairPathsSortedList) {
				debugMes("\t" + p, 15);
			}
			
			
		}
		
		
		//--------------------
		// start assembling
		
		HashMap<List<Integer>, Pair<Integer>> final_transcripts = new HashMap<List<Integer>, Pair<Integer>>();
		
		
		int round = 0;
		while (! pairPathsSortedList.isEmpty()) {
			
			round++;
			debugMes("\n\nPasaFlyUnique, Round: " + round, 10);
			
			ArrayList<PasaVertex> pasaVerticesSortedList = new ArrayList<PasaVertex>();
			for (PairPath pp : pairPathsSortedList) { // already sorted
				int count = pairPathToReadSupport.get(pp);
				pasaVerticesSortedList.add(new PasaVertex(pp, count));
			}

			PasaVertex [] pasaVerticesSortedArr = pasaVerticesSortedList.toArray(new PasaVertex[pasaVerticesSortedList.size()]);
			PairPath[] pairPathsSortedArr = pairPathsSortedList.toArray(new PairPath[pairPathsSortedList.size()]);
			
			
			//----------------------
			// EXAMINE CONTAINMENTS
			//----------------------
			
			// init
			ArrayList<PairPath> pairPathsContainmentsRemoved = new ArrayList<PairPath>(pairPathsSortedList);
			ArrayList<PasaVertex> pasaVerticesContainmentsRemoved = new ArrayList<PasaVertex>(pasaVerticesSortedList);

			
			debugMes("Assigning pairpath containments.", 10);
			List<Integer> containments = assignPasaPairPathContainments(graph, dijkstraDis, pasaVerticesSortedArr); // vertices updated to include containment info.

			debugMes("REMOVING CONTAINMENTS: " + containments, 10);
			for(int i = 0; i < containments.size(); i++)
			{

				pasaVerticesContainmentsRemoved.remove(pasaVerticesSortedArr[containments.get(i)]);
				pairPathsContainmentsRemoved.remove(pairPathsSortedArr[containments.get(i)]);
			}
			
			//------------------------------------------------
			// EXAMINE UNCERTAINTIES THAT BREAK TRANSITIVITY
			//-------------------------------------------

			PairPath [] pairPathsContainmentsRemovedArr = pairPathsContainmentsRemoved.toArray(new PairPath[pairPathsContainmentsRemoved.size()]);
			PasaVertex [] pasaVerticesContainmentsRemovedArr = pasaVerticesContainmentsRemoved.toArray(new PasaVertex[pasaVerticesContainmentsRemoved.size()]);
			
			boolean[][] dag = getPairPathConsistencyDAG(graph, dijkstraDis, pairPathsContainmentsRemovedArr);

			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 10) {
				debugMes("PASA Consistency DAG--------------------------\n", 10);

				System.out.println(boolean_matrix_toString(dag));
			}

			

			ArrayList<PairPath> pairPathsUncertainRemoved = new ArrayList<PairPath>(pairPathsContainmentsRemoved);
			ArrayList<PasaVertex> pasaVerticesUncertainRemoved = new ArrayList<PasaVertex>(pasaVerticesContainmentsRemoved);

		
			debugMes("Identifying uncertain entries that break transitivities.", 10);
			// identify and remove uncertain entries (those that break transitive compatibility relationships)
			ArrayList<Integer> uncertain = getUncertainRequireOverlap(dag, pairPathsContainmentsRemovedArr, graph, dijkstraDis);

			debugMes("Uncertain indices include: " + uncertain, 10);

			
			debugMes("REMOVING UNCERTAINTIES: " + uncertain, 10);

			for(int i = 0; i < uncertain.size(); i++)
			{

				pasaVerticesUncertainRemoved.remove(pasaVerticesContainmentsRemovedArr[uncertain.get(i)]);
				pairPathsUncertainRemoved.remove(pairPathsContainmentsRemovedArr[uncertain.get(i)]);
			}

			PasaVertex[] pasaVerticesUncertainRemovedArr = pasaVerticesUncertainRemoved.toArray(new PasaVertex[pasaVerticesUncertainRemoved.size()]);
			PairPath[] pairPathsUncertainRemovedArr = pairPathsUncertainRemoved.toArray(new PairPath[pairPathsUncertainRemoved.size()]);

			//print pair paths
			debugMes("PAIR PATHS remaining after uncertainties removed ---------------------------------------",10);
			for(int i = 0; i < pairPathsUncertainRemovedArr.length; i++)
			{
				debugMes("PairPathAfterUncertainRemoved "+ i + " " + pairPathsUncertainRemovedArr[i], 10);
			}

			// regenerate the dag now that the uncertain entries are removed.
			dag = getPairPathConsistencyDAG(graph, dijkstraDis, pairPathsUncertainRemovedArr); // already identified containments

			//print dag
			debugMes("DAG after uncertainties removed ---------------------------\n" + boolean_matrix_toString(dag),10);

			// examine neighboring DAG
			for (int i = 0; i < pairPathsUncertainRemovedArr.length-1; i++) {
				if (! dag[i][i+1]) {
					debugMes("NeighborDagCheck: PairPath: [" + i + "] "+ pairPathsUncertainRemovedArr[i] 
							+ "\n\tnot compatible with: [" + (i+1) + "] " + pairPathsUncertainRemovedArr[i+1], 10);
				}
			}
			
			

			//2.2 check transitivity
			if(!checkTransitivityRequireOverlap(dag, pairPathsUncertainRemovedArr, graph, dijkstraDis))
			{
				throw(new RuntimeException("Graph is NOT transitive!"));

			}
			else {
				debugMes("Transitivity of compatibility graph validates.", 10);
			}

			
			debugMes("build_PASA_trellis_left_to_right()", 10);

			build_PASA_trellis_left_to_right(pasaVerticesUncertainRemovedArr, dag, graph, componentReadHash, dijkstraDis, pairPathToReadSupport, tripletMapper, extendedTripletMapper);

			// get highest scoring path:
			debugMes("Identifying highest scoring PASA path.", 10);
			ScoredPath best = null;
			for (int i = 0; i < pasaVerticesUncertainRemovedArr.length; i++) {

				ScoredPath sp = pasaVerticesUncertainRemovedArr[i].get_highest_scoring_fromPath();
				if (best == null || sp.score > best.score) {
					best = sp;
				}
			}


			debugMes("Best score: " + best.score + ", containing entries: " + best.paths, 10);

			List<Integer> best_path_vertex_list = Path.collapse_compatible_pair_paths(best.paths);


			final_transcripts.put(best_path_vertex_list, new Pair<Integer>(1,1));
			debugMes("Reconstructed path is: " + best_path_vertex_list, 10);

			// remove those pairpaths included in the best path
			List<PairPath> toRemove = new ArrayList<PairPath>();
			for (PairPath pp : pairPathsSortedList) {

				if (pp.isCompatibleAndContainedBySinglePath(best_path_vertex_list)) {
					toRemove.add(pp);
					debugMes("compatibly_contained_by_reconstructed_path: " + pp, 10);
				}
				else {
					debugMes("NotCompatibleRetainedForNextRound: " + pp, 10);
				}
			}
			for (PairPath pp : toRemove)
				pairPathsSortedList.remove(pp);

		}
		
		
		return(final_transcripts);
		
	}

	
	
	
	
	
	private static List<PairPath> find_paired_paths_with_greatest_map_support(
			List<ScoredPath> sp_from, List<ScoredPath> sp_to,
			HashMap<PairPath, Integer> finalVertexPositions) {

		
		int best_count = -1;
		List<PairPath> combined = new ArrayList<PairPath>();
		
		for (ScoredPath spA : sp_from) {
		
			HashSet<PairPath> hA = new HashSet<PairPath>();
			for (PairPath h : spA.paths) {
				if (finalVertexPositions.containsKey(h)) {
					hA.add(h);
				}
			}
			
			for (ScoredPath spB : sp_to) {
			
				HashSet<PairPath> hB = new HashSet<PairPath>(hA);
				for (PairPath b : spB.paths) {
					if (finalVertexPositions.containsKey(b)) {
						hB.add(b);
					}
				}
				
				if (hB.size() > best_count) {
					best_count = hB.size();
					combined.clear();
					combined.addAll(spA.paths);
					combined.addAll(spB.paths);
				}
				
				
			}
		}
		
		return(combined);
		
	}


	private static void build_PASA_trellis_left_to_right(
			PasaVertex[] pasaVerticesArr, 
			boolean[][] dag,
			final DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis, 
			Map<PairPath, Integer> pairPathToReadSupport, HashMap<Integer, List<List<Integer>>> tripletMapper, HashMap<Integer, List<List<Integer>>> extendedTripletMapper
			) {
		
		
		
		for  (int i = 1; i < pasaVerticesArr.length; i++)  {
			
			PasaVertex iV = pasaVerticesArr[i];
			
			for (int j = i - 1; j >= 0; j--) {
				
				
				
				PasaVertex iJ = pasaVerticesArr[j];
				
				
				if (! dag[j][i]) {
					// see if too far apart
					if (twoPairPathsAreTooFarAwayInGraph(iV.pp, iJ.pp, graph)) {
						if (FAST_PASA)
							break;
					}
					else {
						continue;  // must conflict
					}
				}
				
				// require that they share a node in common
				if (! iJ.pp.haveAnyNodeInCommon(iV.pp))
					continue;
				
				// see if we can extend paths in iJ to include pairpath represented by iV
				
				final List<ScoredPath> sp_list = iJ.get_fromPaths();
				for (ScoredPath sp : sp_list) {
					// is there sufficient read support for extending this path?
					
					//debugMes("\nnote, sp_list is of size: " + sp_list.size(), 10);
					//debugMes("\nAdding path list to [iV] from [iJ] " + sp.paths, 10);
					
					List<PairPath> extendedList = new ArrayList<PairPath>();
					extendedList.addAll(sp.paths);
					extendedList.add(iV.pp);
					
					if (! violates_triplet_support(tripletMapper, extendedList)) {
						iV.push_fromPaths(new ScoredPath(extendedList, (sp.score + iV.readSupport + iV.num_contained)));
						sp.path_extended = true;
					}
					
				}
				
				
				
			}
			
			
			
		}
		
	}

	
	private static boolean violates_triplet_support(
			HashMap<Integer, List<List<Integer>>> tripletMapper,
			List<PairPath> extendedList) {
		
		
		HashMap<Integer,Integer> prev_node = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> next_node = new HashMap<Integer,Integer>();
		
		for (PairPath pp : extendedList) {
			Integer prev = -1;
			for (List<Integer> path :  pp.get_paths()) {
				for (Integer node : path) {
					if (prev >= 0) {
						prev_node.put(node,  prev);
						next_node.put(prev, node);
					}
					prev = node;
				}
			}

		}
		
		for (Integer center : prev_node.keySet()) {
			if (tripletMapper.containsKey(center) && next_node.containsKey(center)) {
				Integer left = prev_node.get(center);
				Integer right = next_node.get(center);
				
				List<Integer> curr_triplet = new ArrayList<Integer>();
				curr_triplet.add(left);
				curr_triplet.add(center);
				curr_triplet.add(right);
				
				List<List<Integer>> triplets = tripletMapper.get(center);
				if (! tripletSupported(triplets, curr_triplet)) {
					debugMes("PASA TRIPLET CHECK WARNING: triplet: " + curr_triplet + " violates available triplets: " + triplets + " and so path list is not valid: " + extendedList, 15);
					return(true); // yes, violates
				}
			}
		}
		
		
		return false; // no violation found.
	}


	private static void build_PASA_trellis_right_to_left (
			PasaVertex[] pasaVerticesArr, 
			boolean[][] dag,
			final DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<Integer, HashMap<PairPath, Integer>> componentReadHash,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis, 
			Map<PairPath, Integer> pairPathToReadSupport, HashMap<Integer, List<List<Integer>>> tripletMapper
			) {
		
	
		
		for  (int i = pasaVerticesArr.length-2; i >= 0; i--)  {
			
			PasaVertex iV = pasaVerticesArr[i];
			
			for (int j = i +1; j < pasaVerticesArr.length; j++) {
				
				
				PasaVertex iJ = pasaVerticesArr[j];
				
				if (! dag[i][j]) {
					if (twoPairPathsAreTooFarAwayInGraph(iV.pp, iJ.pp, graph)) {
						if (FAST_PASA)
							break;
					}
					else {
						continue;  // must conflict
					}
					
				}
				
				
				
				// require that they share a node in common
				if (! iJ.pp.haveAnyNodeInCommon(iV.pp))
					continue;
				
				// see if we can extend paths in iJ to include pairpath represented by iV
				
				final List<ScoredPath> sp_list = iJ.get_toPaths();
				for (ScoredPath sp : sp_list) {
					// is there sufficient read support for extending this path?
					//   * implement later on * //
					
					List<PairPath> extendedList = new ArrayList<PairPath>();
					extendedList.addAll(sp.paths);
					extendedList.add(iV.pp);
					
					if (! violates_triplet_support(tripletMapper, extendedList)) {
						iV.push_toPaths(new ScoredPath(extendedList, (sp.score + iV.readSupport + iV.num_contained)));
						
					}
					
				}
				
				
				
			}
			
			
			
		}
		
		
		
	}

	
	
	
	
	private static String boolean_matrix_toString (boolean [][] dag) {
		
		//print dag
		
		if (dag.length > 200) {
			debugMes("dag matrix too large to print in a useful way.", 10);
			return("");
		}
		
		String dag_text = "";
		for(int i = 0; i < dag.length; i++)
		{
			
			for(int j = 0; j < dag.length; j++)
			{
				
				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15)
					System.err.print("\r[" + i + "," + j + "]   " );

				dag_text += ((dag[i][j]) ? 1: 0)  + " ";
			}
			dag_text += "\n";
		}
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 10)
			System.err.println();
		
		return(dag_text);
	}
	

	private static List<Integer> assignPasaPairPathContainments(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis,
			PasaVertex[] pasaVerticesArr
			
			) {
		
		
		
		HashMap<Integer,Boolean> containments = new HashMap<Integer,Boolean>();
		
		for (int i = 0; i < pasaVerticesArr.length; i++) {
			
			PasaVertex iV = pasaVerticesArr[i];
			
			for (int j = 0; j < pasaVerticesArr.length; j++) {
				
				if (i==j)
					continue;
					
				PasaVertex iJ = pasaVerticesArr[j];
		
				if (iJ.pp.haveAnyNodeInCommon(iV.pp) && (iV.pp.isCompatibleAndContainedByPairPath(iJ.pp, graph, dijkstraDis))) {
						
					iJ.num_contained += iV.readSupport;
					containments.put(i, true);
					debugMes("Containment: " + iV.pp + " is contained by: " + iJ.pp, 10);
				}
				
			}
			
			
		}
		
		List<Integer> containment_list = new ArrayList<Integer>(containments.keySet());
		
		return(containment_list);
		
	}


	private static HashMap<List<Integer>, HashMap<PairPath, Integer>> assignCompatibleReadsToPaths(
			HashMap<List<Integer>, Pair<Integer>> finalPaths_all,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash) {
		
		debugMes("\n\n## assignCompatibleReadsToPaths()", 20);
		
		HashMap<List<Integer>, HashMap<PairPath, Integer>> pathToContainedReads = new HashMap<List<Integer>, HashMap<PairPath, Integer>>();
		
		for (List<Integer> path : finalPaths_all.keySet()) {
			
			for (HashMap<PairPath,Integer> read_map : combinedReadHash.values()) {
				for (PairPath p : read_map.keySet()) {
					
					if (p.isCompatibleAndContainedBySinglePath(path)) {
						if (! pathToContainedReads.containsKey(path)) {
							pathToContainedReads.put(path, new HashMap<PairPath, Integer>());
							
						}
						
						debugMes("assignCompatibleReadsToPaths: " + p + " is compatible with " + path, 20);
						
						pathToContainedReads.get(path).put(p, read_map.get(p));
					}
					else {
						debugMes("assignCompatibleReadsToPaths: " + p + " is NOT compatible with " + path, 20);
					}
					
				}
			}
			
		}
		
		return (pathToContainedReads);
	}


	private static boolean containsNull(PairPath pp1)
	{
		if(pp1.getFirstID() == null)
			return true;
		else
			return false;
	}

	/**
	 * Given two pair paths, determines consistency/compatibility of two pair paths for the partial order.
	 * @param pp1
	 * @param pp2
	 * @param graph
	 * @param dijkstraDis
	 * @return 0 if they are not consistent, 1 if pp1 comes before pp2, -1 if pp2 comes before pp1
	 */
	private static int isConsistent(PairPath pp1, PairPath pp2, DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis)
	{
		
		pp1 = pp1.trimSinkNodes();
		pp2 = pp2.trimSinkNodes();
		
		debugMes("isConsistent? " + pp1 + pp2, 15);
		
		if (pp1.equals(pp2)) { return (0); }
		
		
		// If have nodes in common, require compatibility in overlapping path parts:
		
		if(pp1.haveAnyNodeInCommon(pp2))
		{
			debugMes("\tHave nodes in common.", 15);
			
			if(!pp1.isCompatible(pp2))
			{
				debugMes("\tNot compatible.", 15);
				return 0;
			}
		}
		
		// iterate through every node in pp1 and check if it's consistent
		List<Integer> path1 = pp1.getPath1();
	
		// see that pp1:path1 nodes are consistent with pp2
		Iterator itr = path1.iterator();
		for(int i = 0; i < path1.size(); i++)
		{
			Integer n = (Integer)(itr.next());
			
			if(!(readIsConsistentWithNode(pp2, n, graph, dijkstraDis))) {
				debugMes("\tpp2: " + pp2 + " is not consistent with node: " + n, 15);
				return 0;
			}
		}
		
			
		if(pp1.hasSecondPath())
		{

			// see if pp1:path2 are consistent with pp2
			List<Integer> path2 = pp1.getPath2();
			itr = path2.iterator();
			for(int i = 0; i < path2.size(); i++)
			{
				Integer n = (Integer)(itr.next());
				if(!(readIsConsistentWithNode(pp2, n, graph, dijkstraDis))) {
					debugMes("\tpp2: " + pp2 + " second path is not consistent with node: " + n, 15);
					return 0;
				}
			}


		}
		

		// require that v1 comes before v2 in the partial order
		SeqVertex v1 = getSeqVertex(graph, pp1.getFirstID());
		SeqVertex v2 = getSeqVertex(graph, pp2.getFirstID());
		
		if (v1.equals(v2)) {
			// check their last nodes.
			SeqVertex lv1 = getSeqVertex(graph, pp1.getLastID());
			SeqVertex lv2 = getSeqVertex(graph, pp2.getLastID());
			
			if (lv1.equals(lv2)) {
				// must have same first and same last node, so potential differences in-between but otherwise compatible afaict
				int pp1_hashcode = pp1.hashCode();
				int pp2_hashcode = pp2.hashCode();
				int consistent = (pp1.hashCode() < pp2.hashCode()) ? 1 : 0;
				
				debugMes("\tfirst vertex node: " + pp1.getFirstID() + " and last node " + pp1.getLastID() + " are equal, so defining consistency based on hashcode comparison.", 15);
				return(consistent); // just use consistent ordering to define proper DAG connectability
			}
			else {
				
				// first node equivalent, last node not equivalent
				
				int ancestral = SeqVertex.isAncestral(lv1, lv2, dijkstraDis);
				debugMes("\tfirst nodes same: " + pp1.getFirstID() + ", but last nodes are different: " + 
						pp1.getLastID() + " vs. " + pp2.getLastID() + ", and SeqVertex.isAncestral = " + ancestral, 15);
				return( (ancestral>0) ? 1:0);
			}
			
		}
		else {
			
			// first node not equivalent.
			
			int ancestral = SeqVertex.isAncestral(v1,v2,dijkstraDis);
			
			debugMes("\tpairpaths are compatible, examining relative orientation of first vertices: " 
					+ v1.getID() + " vs. " + v2.getID() + ", ancestral = " + ancestral, 15);
			return( (ancestral > 0) ? 1 : 0);
		}
		
		
		
		
		
	}
	
	
	private static boolean isOverlappingAndDirectionallyConsistent(PairPath pp1, PairPath pp2, DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis)
	{
		
		pp1 = pp1.trimSinkNodes();
		pp2 = pp2.trimSinkNodes();
		
		debugMes("isOverlappingAndDirectionallyConsistent? " + pp1 + pp2, 15);
		
		if (pp1.equals(pp2)) { return (true); }
		
		
		// If have nodes in common, require compatibility in overlapping path parts:
		
		if(pp1.haveAnyNodeInCommon(pp2))
		{
			debugMes("\tHave nodes in common.", 15);
			
			if(!pp1.isCompatible(pp2))
			{
				debugMes("\tNot compatible.", 15);
				return false;
			}
			
			// DO allow containments to be compatible
			if (pp2.isCompatibleAndContainedByPairPath(pp1, graph, dijkstraDis)) {
				debugMes("\tpp2 isCompatibleAndContainedBy pp1, setting true (containments removed later on).", 15);
				return(true);
			}
		}
		else {
			debugMes("\tNo node overlap, so not compatible.", 15);
			return(false);
		}
		
		// iterate through every node in pp1 and check if it's consistent
		List<Integer> path1 = pp1.getPath1();
	
		// see that pp1:path1 nodes are consistent with pp2
		Iterator<Integer> itr = path1.iterator();
		for(int i = 0; i < path1.size(); i++)
		{
			Integer n = itr.next();
			
			if(!(readIsConsistentWithNode(pp2, n, graph, dijkstraDis))) {
				debugMes("\tpp2: " + pp2 + " is not consistent with node: " + n, 15);
				return false;
			}
		}
		
			
		if(pp1.hasSecondPath())
		{

			// see if pp1:path2 are consistent with pp2
			List<Integer> path2 = pp1.getPath2();
			itr = path2.iterator();
			for(int i = 0; i < path2.size(); i++)
			{
				Integer n = (Integer)(itr.next());
				if(!(readIsConsistentWithNode(pp2, n, graph, dijkstraDis))) {
					debugMes("\tpp2: " + pp2 + " second path is not consistent with node: " + n, 15);
					return false;
				}
			}


		}

		// require that v1 comes before v2 in the partial order
		SeqVertex v1 = getSeqVertex(graph, pp1.getFirstID());
		SeqVertex v2 = getSeqVertex(graph, pp2.getFirstID());

		if (v1.equals(v2)) {
			// check their last nodes.
			SeqVertex lv1 = getSeqVertex(graph, pp1.getLastID());
			SeqVertex lv2 = getSeqVertex(graph, pp2.getLastID());

			if (lv1.equals(lv2)) {
				// must have same first and same last node, so potential differences in-between but otherwise compatible afaict
				return(true);
			}
			else {

				// first node equivalent, last node not equivalent

				int ancestral = SeqVertex.isAncestral(lv1, lv2, dijkstraDis);
				debugMes("\tfirst nodes same: " + pp1.getFirstID() + ", but last nodes are different: " + 
						pp1.getLastID() + " vs. " + pp2.getLastID() + ", and SeqVertex.isAncestral = " + ancestral, 15);
				return(ancestral>0);
			}

		}
		else {

			// first node not equivalent.

			int ancestral = SeqVertex.isAncestral(v1,v2,dijkstraDis);

			debugMes("\tpairpaths are compatible, examining relative orientation of first vertices: " 
					+ v1.getID() + " vs. " + v2.getID() + ", ancestral = " + ancestral, 15);
			return(ancestral > 0);
		}
		
		
	}
	
	
	private static boolean checkTransitivity(int[][] adj, PairPath[] pairPathArr, HashMap<PairPath,Integer> pairPathToIntVal)
	{
	
		// Examine all triplets (A,B,C) and ensure that A->B and B->C implies A->C
		
		
		for(int i = 0; i < adj.length; i++)
		{

			for(int j = 0; j < adj.length; j++)
			{
					
				for(int k = 0; k < adj.length; k++)
				{
					
					if (adj[i][j] == 1 && adj[j][k] == 1 && adj[i][k] == 0) {
						
						Integer orig_i = pairPathToIntVal.get(pairPathArr[i]);
						Integer orig_j = pairPathToIntVal.get(pairPathArr[j]);
						Integer orig_k = pairPathToIntVal.get(pairPathArr[k]);
						

							debugMes("UNCERTAINTY DETECTED AFTER SUPPOSEDLY HAVING REMOVED THEM [" + orig_i + "," + orig_j + "," + orig_k + "] :\n" +
									i + " " + pairPathArr[i] + " is consistent with " + j + " " + pairPathArr[j] + "\n" +
									j + " " + pairPathArr[j] + " is consistent with " + k + " " + pairPathArr[k] + "\n" +
									i + " " + pairPathArr[i] + " is NOT consistent with " + k + " " + pairPathArr[k] + "\n", 10);
						
						return false;
				
					}
					
				}
			}
		}
		return true;
	}
	
	private static boolean checkTransitivity(boolean[][] adj, PairPath[] pairPathArr, HashMap<PairPath,Integer> pairPathToIntVal)
	{
	
		// Examine all triplets (A,B,C) and ensure that A->B and B->C implies A->C
		
		
		for(int i = 0; i < adj.length; i++)
		{

			for(int j = 0; j < adj.length; j++)
			{
					
				for(int k = 0; k < adj.length; k++)
				{
					
					if (adj[i][j] && adj[j][k] && ! adj[i][k]) {
						
						Integer orig_i = pairPathToIntVal.get(pairPathArr[i]);
						Integer orig_j = pairPathToIntVal.get(pairPathArr[j]);
						Integer orig_k = pairPathToIntVal.get(pairPathArr[k]);
						

							debugMes("UNCERTAINTY DETECTED AFTER SUPPOSEDLY HAVING REMOVED THEM [" + orig_i + "," + orig_j + "," + orig_k + "] :\n" +
									i + " " + pairPathArr[i] + " is consistent with " + j + " " + pairPathArr[j] + "\n" +
									j + " " + pairPathArr[j] + " is consistent with " + k + " " + pairPathArr[k] + "\n" +
									i + " " + pairPathArr[i] + " is NOT consistent with " + k + " " + pairPathArr[k] + "\n", 10);
						
						return false;
				
					}
					
				}
			}
		}
		return true;
	}
	
	private static boolean checkTransitivityRequireOverlap(
			boolean[][] adj, 
			PairPath[] pairPathArr, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis)
	{
	
		// Examine all triplets (A,B,C) and ensure that A->B and B->C implies A->C
		
		
		for(int i = 0; i < adj.length-2; i++)
		{

			for(int j = i+1; j < adj.length-1; j++)
			{
				
				// see if j is too far away from i, can then go to next i

				if (! adj[i][j]) {
					if (twoPairPathsAreTooFarAwayInGraph(pairPathArr[i], pairPathArr[j], graph))	
						if (FAST_PASA)
							break;
				}
				
				
					
				for(int k = j+1; k < adj.length; k++)
				{
					
					
					// see if k is too far away from j, can then go to next j
					if (! adj[j][k]) {
						if (twoPairPathsAreTooFarAwayInGraph(pairPathArr[j], pairPathArr[k], graph))
							if (FAST_PASA)
								break;
					}
					
					
					
					if (adj[i][j] && adj[j][k] && pairPathArr[i].haveAnyNodeInCommon(pairPathArr[k]) && ! adj[i][k]) {
						
						

							debugMes("UNCERTAINTY DETECTED AFTER SUPPOSEDLY HAVING REMOVED THEM [" + i + "," + j + "," + k + "] :\n" +
									i + " " + pairPathArr[i] + " is consistent with " + j + " " + pairPathArr[j] + "\n" +
									j + " " + pairPathArr[j] + " is consistent with " + k + " " + pairPathArr[k] + "\n" +
									i + " " + pairPathArr[i] + " is NOT consistent with " + k + " " + pairPathArr[k] + "\n", 10);
						
						return false;
				
					}
					
				}
			}
		}
		return true;
	}
	
	private static ArrayList<Integer> getUncertain(int[][] adj, PairPath[] pairPathArr)
	{
		
		// DAG:   i -> j -> k
		
	
		ArrayList<Integer> toRemove = new ArrayList<Integer>();
		
		for(int i = 0; i < adj.length; i++)
		{
			if (toRemove.contains(i)) { continue; }
			
		
			// move j to previous shared first node with i in the list:
			int j = i;
			if (j > 0 ) {
				while (j > 0 && pairPathArr[i].getFirstID().equals(pairPathArr[j-1].getFirstID())) {
					j--;
				}
			}
			
			
			
			for(; j < adj.length; j++)
			{
				
				
				if(toRemove.contains(j))
					continue;
		
				// start k at first node sharing the same first ID as j
				int k = j;
				if (k > 0) {
					while (k>0 && pairPathArr[j].getFirstID().equals(pairPathArr[k-1].getFirstID())) {
						k--;
					}
				}
				
				
				for(; k < adj.length; k++)
				{
					if(toRemove.contains(k))
						continue;
					
					debugMes("CHECKING TRANSITIVITY [" + i + "," + j + "," + k + "] "
							+ "= [" + adj[i][j] + "," + adj[j][k] + "," +  adj[i][k] + "]", 15);
						
					
					debugMes("MORE VERBOSE CHECKING TRANSITIVITY:] " +
							" { " + i + " " + pairPathArr[i] + " results(" + adj[i][j] + ") " + j + " " + pairPathArr[j] + " }" +
							" { " + j + " " + pairPathArr[j] + " results("+ adj[j][k] + ") " + k + " " + pairPathArr[k] + " } " +
							" { " + i + " " + pairPathArr[i] + " results(" + adj[i][k] + ") " + k + " " + pairPathArr[k] + " } ", 18);
					
					
					if (adj[i][j] == 1 && adj[j][k] == 1)
					{
					
						
						if (adj[i][k] == 0) {
						
							toRemove.add(j); // central node breaks transitivity.  Remove it. 
					

							debugMes("UNCERTAINTY DETECTED:, removing: " +  pairPathArr[j] + "\n" +
									i + " " + pairPathArr[i] + " is consistent with " + j + " " + pairPathArr[j] + "\n" +
									j + " " + pairPathArr[j] + " is consistent with " + k + " " + pairPathArr[k] + "\n" +
									i + " " + pairPathArr[i] + " is NOT consistent with " + k + " " + pairPathArr[k] + "\n", 10);

							break; // go to next j
						}
					}
				}
			}
		}
		
		return toRemove;
	}
	
	private static ArrayList<Integer> getUncertain(boolean[][] adj, PairPath[] pairPathArr)
	{
		
		// DAG:   i -> j -> k
		
		
		int removeFlag = 0;
		ArrayList<Integer> toRemove = new ArrayList<Integer>();
		
		for(int i = 0; i < adj.length; i++)
		{
			if (toRemove.contains(i)) { continue; }
			
		
			// move j to previous shared first node with i in the list:
			int j = i;
			if (j > 0 ) {
				while (j > 0 && pairPathArr[i].getFirstID().equals(pairPathArr[j-1].getFirstID())) {
					j--;
				}
			}
			
			
			
			for(; j < adj.length; j++)
			{
				
				
				if(toRemove.contains(j))
					continue;
		
				// start k at first node sharing the same first ID as j
				int k = j;
				if (k > 0) {
					while (k>0 && pairPathArr[j].getFirstID().equals(pairPathArr[k-1].getFirstID())) {
						k--;
					}
				}
				
				
				for(; k < adj.length; k++)
				{
					if(toRemove.contains(k))
						continue;
					
					if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15)
						System.err.print("\r[" + i + "," + j + "," + k + "]  ");
					
					debugMes("CHECKING TRANSITIVITY [" + i + "," + j + "," + k + "] "
							+ "= [" + adj[i][j] + "," + adj[j][k] + "," +  adj[i][k] + "]", 15);
						
					
					debugMes("MORE VERBOSE CHECKING TRANSITIVITY:] " +
							" { " + i + " " + pairPathArr[i] + " results(" + adj[i][j] + ") " + j + " " + pairPathArr[j] + " }" +
							" { " + j + " " + pairPathArr[j] + " results("+ adj[j][k] + ") " + k + " " + pairPathArr[k] + " } " +
							" { " + i + " " + pairPathArr[i] + " results(" + adj[i][k] + ") " + k + " " + pairPathArr[k] + " } ", 18);
					
					
					if (adj[i][j] == true && adj[j][k] == true)
					{
					
						
						if (adj[i][k] == false) {
						
							toRemove.add(j); // central node breaks transitivity.  Remove it. 
					

							debugMes("UNCERTAINTY DETECTED:, removing: " +  pairPathArr[j] + "\n" +
									i + " " + pairPathArr[i] + " is consistent with " + j + " " + pairPathArr[j] + "\n" +
									j + " " + pairPathArr[j] + " is consistent with " + k + " " + pairPathArr[k] + "\n" +
									i + " " + pairPathArr[i] + " is NOT consistent with " + k + " " + pairPathArr[k] + "\n", 10);

							break; // go to next j
						}
					}
				}
			}
		}
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 12) 
			System.err.println();
		
		return toRemove;
	}
	
	
	private static ArrayList<Integer> getUncertainRequireOverlap(
			boolean[][] adj, 
			PairPath[] pairPathArr,
			final DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis)
	{
		
		// DAG:   i -> j -> k
		
		
		int debug_node = 301;
		boolean local_debug = false;
		
		ArrayList<Integer> toRemove = new ArrayList<Integer>();
		
		for(int i = 0; i < adj.length-2; i++)
		{
			
			if (local_debug && i == debug_node)
				System.err.print("\rUncertaintyCheck: [" + i + "," + "?" + "," + "?" + "]\n" );
			
						
			if (toRemove.contains(i)) { 
				if (local_debug && i == debug_node)
					System.err.print("\rUncertaintyCheck: [" + i + "," + "?" + "," + "?" + "]  Iprev-Captured\n");
				
				continue; 
			}
			
		

			boolean tooFar = false;
			
			for(int j = i + 1; j < adj.length-1; j++)
			{
				
				if (local_debug && i == debug_node)
					System.err.print("\rUncertaintyCheck: [" + i + "," + j + "," + "?" + "]\n");
				
				
				if(toRemove.contains(j)) {
				
					if (local_debug && i == debug_node)
						System.err.print("\rUncertaintyCheck: [" + i + "," + j + "," + "?" + "]  Jprev-Captured\n");
					
					continue;
				}
				
				// see if j is too far away from i, can then go to next i

				if (! adj[i][j]) {
					if (twoPairPathsAreTooFarAwayInGraph(pairPathArr[i], pairPathArr[j], graph)) {	
						if (FAST_PASA)
							break;
					}
				}
				
				if (twoPairPathsAreTooFarAwayInGraph(pairPathArr[i], pairPathArr[j], graph))	
					tooFar = true;
				else if (tooFar)
						debugMes("CHANGED from too far to within distance again: I:" + pairPathArr[i] + " J:" + pairPathArr[j], 10);


				
				boolean tooFar2 = false;
				
				for(int k = j + 1; k < adj.length; k++)
				{
					
					if (local_debug && i == debug_node)
						System.err.print("\rUncertaintyCheck: [" + i + "," + j + "," + k + "] " +
					    "[" + adj[i][j] + "," + adj[j][k] + "," + adj[i][k] + "]\n" );
					
					
					
					if(toRemove.contains(k)) {
						if (local_debug && i == debug_node)
							System.err.print("\rUncertaintyCheck: [" + i + "," + j + "," + k + "]  Kprev-Captured\n");
						
						continue;
					}

					// see if k is too far away from j, can then go to next j
					if (! adj[j][k]) {
						if (twoPairPathsAreTooFarAwayInGraph(pairPathArr[j], pairPathArr[k], graph)) {
							if (FAST_PASA)
								break;
						}
					}

					if (twoPairPathsAreTooFarAwayInGraph(pairPathArr[j], pairPathArr[k], graph))
						tooFar2 = true;
					else if (tooFar2)
							debugMes("CHANGED from too far to within distance again: I:" + pairPathArr[i] + " J:" + pairPathArr[j], 10);

					
					
					if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
						System.err.print("\r[" + i + "," + j + "," + k + "]  ");
					

						debugMes("CHECKING TRANSITIVITY [" + i + "," + j + "," + k + "] "
								+ "= [" + adj[i][j] + "," + adj[j][k] + "," +  adj[i][k] + "]", 15);


						debugMes("MORE VERBOSE CHECKING TRANSITIVITY:] " +
								" { " + i + " " + pairPathArr[i] + " results(" + adj[i][j] + ") " + j + " " + pairPathArr[j] + " }" +
								" { " + j + " " + pairPathArr[j] + " results("+ adj[j][k] + ") " + k + " " + pairPathArr[k] + " } " +
								" { " + i + " " + pairPathArr[i] + " results(" + adj[i][k] + ") " + k + " " + pairPathArr[k] + " } ", 18);

					}

					if (adj[i][j] == true && adj[j][k] == true)
					{


						if (pairPathArr[i].haveAnyNodeInCommon(pairPathArr[k]) && adj[i][k] == false) {

							toRemove.add(j); // central node breaks transitivity.  Remove it. 

							

							debugMes("UNCERTAINTY DETECTED:, removing: " +  pairPathArr[j] + "\n" +
									i + " " + pairPathArr[i] + " is consistent with " + j + " " + pairPathArr[j] + "\n" +
									j + " " + pairPathArr[j] + " is consistent with " + k + " " + pairPathArr[k] + "\n" +
									i + " " + pairPathArr[i] + " is NOT consistent with " + k + " " + pairPathArr[k] + "\n", 10);

							
							if (pairPathArr[j].isCompatibleAndContainedByPairPath(pairPathArr[i]) || pairPathArr[j].isCompatibleAndContainedByPairPath(pairPathArr[k])) {
								throw new RuntimeException("ERROR, containment encountered where containments should have been removed.");
							}
							
							
							if (local_debug && i == debug_node)
								System.err.print("\rUncertaintyCheck: [" + i + "," + j + "," + k + "]  J-Captured\n");
							
							break; // no more need to analyze k in this i-j-k series.
						}
						else {
							if (local_debug && i == debug_node)
								System.err.print("\rUncertaintyCheck: [" + i + "," + j + "," + k + "]  OK\n");
						}
					}
				}
			}
			
	
		}
		
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 12) 
			System.err.println();
		
		
		return toRemove;
	}
	
	
	
	private static int[][] getPairPathDAG(DirectedSparseGraph<SeqVertex, SimpleEdge> graph, DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis,
			PairPath[] pairPathArr)
	{
		
		debugMes("getPairPathDAG:", 10);
		
		int[][] dag = new int[pairPathArr.length][pairPathArr.length];
		
		for (int[] row : dag)
		    Arrays.fill(row, 0); // init to no connection.
		
	
		for (int i = 0; i < pairPathArr.length; i++)
		{
			if (pairPathArr[i].isEmpty())
				continue;
			
			//start comparisons to j where j starts at least at the same position as path i.
			int j = i;
			
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15 && j>0) {
				Integer i_first_id = pairPathArr[i].getFirstID();
				Integer jm1_first_id = pairPathArr[j-1].getFirstID();
			
				debugMes("-comparing first IDs for :[" + i +"," + j + "-1]: " + i_first_id + " to " + jm1_first_id, 15);
			}
			
			while (j > 0 && pairPathArr[i].getFirstID().equals(pairPathArr[j-1].getFirstID())) {
				
				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
					Integer i_first_id = pairPathArr[i].getFirstID();
					Integer jm1_first_id = pairPathArr[j-1].getFirstID();
					debugMes("-comparing first IDs for :[" + i +"," + j + "-1]: " + i_first_id + " to " + jm1_first_id, 15);
				}
				j--;
			}
			
			for (; j < pairPathArr.length; j++)
			{
				
				int value;
				
				if (i == j) 
					value = 1;  // make compatible for now.  Remove self-compatibility before maximal matching
				else {
					value = isConsistent(pairPathArr[i], pairPathArr[j], graph, dijkstraDis);
					
					// nope, below doesn't work - maximal matching does require the full set of compatibilities.
					//boolean val = isOverlappingAndDirectionallyConsistent(pairPathArr[i], pairPathArr[j], graph, dijkstraDis);
					//value = (val) ? 1 : 0;
				}
				
				dag[i][j] = value;
				
				debugMes("Comparing node " + i +" " + pairPathArr[i] + " with node " +  j +" " + pairPathArr[j] + "Result: " + dag[i][j],15);
				
				debugMes("DAG[" + i + "," + j + "]=" + dag[i][j], 15);
				
				if (j < i && value == 1 && ! pairPathArr[i].getFirstID().equals(pairPathArr[j].getFirstID())) {
					debugMes("\tWARNING: ConsistencyListUnordered: [" + i + "," + j + "] " + pairPathArr[i] +  pairPathArr[j], 10); // perhaps should be more serious - throw exception?
					//throw(new RuntimeException("ERROR: ConsistencyListUnordered: [" + i + "," + j + "] " + pairPathArr[i] +  pairPathArr[j]));
				}
			
			}
		}
		return dag;
	}
	
	
	private static boolean[][] getPairPathConsistencyDAG(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis,
			PairPath[] pairPathArr)
	{
		
		debugMes("getPairPathDAG:", 10);
		
		boolean[][] dag = new boolean[pairPathArr.length][pairPathArr.length];
		
		for (boolean[] row : dag)
		    Arrays.fill(row, false); // init to no connection.
		
		
		// i -> j
	
		for (int i = 0; i < pairPathArr.length-1; i++)
		{
			
			PairPath pp_i = pairPathArr[i];
			
			
			boolean tooFar = false;
			 
			for (int j = i + 1; j < pairPathArr.length; j++)
			{
				
				PairPath pp_j = pairPathArr[j];
				
				
				 boolean compatible = isOverlappingAndDirectionallyConsistent(pp_i, pp_j, graph, dijkstraDis);

				 dag[i][j] =compatible;
				 
				
				 if (twoPairPathsAreTooFarAwayInGraph(pp_i, pp_j, graph) && compatible) {
					 debugMes("HOW CAN THESE BE TOO FAR AWAY AND STILL COMPATIBLE? " + pp_i + " vs. " + pp_j, 10);
					 debugMes(report_node_depths(pp_i, graph), 10);
					 debugMes(report_node_depths(pp_j, graph), 10);
					 
					 
				 }
				 if (! compatible) {
					
					 if (twoPairPathsAreTooFarAwayInGraph(pp_i, pp_j, graph)) {
						 if (FAST_PASA)
							 break; 
						 
					 }
				 }
				 
				 if (twoPairPathsAreTooFarAwayInGraph(pp_i, pp_j, graph)) {
					 tooFar = true;
				 }
				 else if (tooFar)
					 debugMes("NOT_TOO_FAR_AFTER_ALL: [" + i + "," + j + "]", 10);
				 
				
				debugMes("Comparing node " + i +" " + pp_i + " with node " +  j +" " + pp_j + "Result: " + compatible,15);
				
				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15)
					System.err.print("\rDAG[" + i + "," + j + "]=" + dag[i][j]);
				
			}
		}
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 12)
			System.err.println();
		
		return dag;
	}
	
	
	private static String report_node_depths(PairPath pp,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
	
		int path_counter = 0;
		String node_depth_text = "";
		for (List<Integer> path : pp.get_paths()) {
			 path_counter++;
			 node_depth_text += "path pt" + path_counter + ": " + path + "\n";
			 for ( Integer node_id : path) {
				 SeqVertex v = getSeqVertex(graph, node_id);
			 
				 node_depth_text += "\tnode: " + node_id + " depth: " + v._node_depth + "\n";
			 }
		}
		
		
		return node_depth_text;
	}
	

	private static ArrayList<Integer> extendChain(ArrayList<Integer> extractedVerticesIDs, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, Map<PairPath, Integer> pairPathToReadSupp, DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis)
	{
		debugMes("Extending Chain", 10);
		ArrayList<Integer> extractedVerticesIDExtended = new ArrayList<Integer>();
		extractedVerticesIDExtended.addAll(extractedVerticesIDs);
		
		Integer firstID = extractedVerticesIDs.get(0);
		Integer lastID = extractedVerticesIDs.get(extractedVerticesIDs.size() - 1);
		
	
		
		// extend from left:
		
		Integer lastIDVisited = firstID;

		
		boolean canExtend = true;
		if(lastIDVisited < 0) // if first ID indicates a source vertex, cannot extend further
			canExtend = false;
		
		// Greedily extend from first ID
		while(canExtend == true)
		{
			int best = 0;
			List<Integer> best_extension_path = null;
			
			canExtend = false;
			firstID = extractedVerticesIDExtended.get(0);
			
			lastIDVisited = firstID;
			
			
			if(lastIDVisited < 0) // if last ID indicates a sink vertex, cannot extend further
				break;
			
			for(PairPath p : pairPathToReadSupp.keySet())
			{
				if(p.isEmpty())
					continue;
				
				int support = pairPathToReadSupp.get(p);
				
				p = p.trimSinkNodes();
				
				if (! p.isCompatible(extractedVerticesIDExtended)) { continue; }
				
				if (! p.containsID(lastIDVisited)) { continue; }
				
				List<Integer> extensionPath;
				if (p.getPath1().contains(lastIDVisited))
					extensionPath = p.getPath1();
				else if (p.getPath2().contains(lastIDVisited))
					extensionPath = p.getPath2();
				else
					throw (new RuntimeException("error, pairpath " + p + " is missing required id " + lastIDVisited));
				
				if(extensionPath.get(0).equals(lastIDVisited)) 
					continue;  // nothign to extend left with since left-most node is same as current path trying to extend.
				
				
				// see if encountering nodes already in the list, indicative of a loop
				Integer loc_in_list = extensionPath.indexOf(lastIDVisited);
				List<Integer> extension_nodes = extensionPath.subList(0, loc_in_list);
				if (PairPath.haveAnyNodeInCommon(extension_nodes, extractedVerticesIDExtended))
					continue;
				
				
				// must have a candidate for extension.
				

				if(support > best)
				{
					best = support;
					best_extension_path = extensionPath;
					canExtend = true;
				}
				
			} // end of testing for extension
			
			if(canExtend == false)
				break;
		
			
			debugMes("Left-extension of : " + extractedVerticesIDExtended + " with " + best_extension_path, 10);
			
			Integer loc_in_list = best_extension_path.indexOf(lastIDVisited);
			
			extractedVerticesIDExtended.addAll(0,  best_extension_path.subList(0, loc_in_list));
		
		}
		//System.exit(0);
		
		// ----------------
		// extend to right
		// ----------------
		
		lastIDVisited = lastID;
		canExtend = true;
		if(lastIDVisited < 0) // if last ID indicates a sink vertex, cannot extend further
			canExtend = false;
		// Greedily extend from last ID
		while(canExtend == true)
		{
			int best = 0;
			List<Integer> best_extension_path = null;
			
			canExtend = false;
			lastID = extractedVerticesIDExtended.get(extractedVerticesIDExtended.size()-1);
			
			lastIDVisited = lastID;
			
			if(lastIDVisited < 0) // if last ID indicates a sink vertex, cannot extend further
				break;
			
			for(PairPath p : pairPathToReadSupp.keySet())
			{
				if(p.isEmpty())
					continue;
				
				int support = pairPathToReadSupp.get(p);
				
				p = p.trimSinkNodes();
				
				if (! p.isCompatible(extractedVerticesIDExtended)) { continue; }
				
				if (! p.containsID(lastIDVisited)) { continue; }
				
				List<Integer> extensionPath;
				if (p.getPath1().contains(lastIDVisited))
					extensionPath = p.getPath1();
				else if (p.getPath2().contains(lastIDVisited))
					extensionPath = p.getPath2();
				else
					throw (new RuntimeException("error, pairpath " + p + " is missing required id " + lastIDVisited));
				
				if(extensionPath.get(extensionPath.size()-1).equals(lastIDVisited))
					continue;
				
				
				// see if encountering nodes already in the list, indicative of a loop
				Integer loc_in_list = extensionPath.indexOf(lastIDVisited);
				List<Integer> extension_nodes = extensionPath.subList(loc_in_list+1, extensionPath.size());
				if (PairPath.haveAnyNodeInCommon(extension_nodes, extractedVerticesIDExtended))
					continue;
				
				
				
				// must have a candidate for extension.
				

				if(support > best)
				{
					best = support;
					best_extension_path = extensionPath;
					canExtend = true;
				}
				
			} // end of testing for extension
			
			if(canExtend == false)
				break;
			
			
			debugMes("Right-extension of : " + extractedVerticesIDExtended + " with " + best_extension_path, 10);
			
			Integer loc_in_list = best_extension_path.indexOf(lastIDVisited);
			
			debugMes("\tloc of " + lastIDVisited + " in best extension list: " + best_extension_path + " is " + loc_in_list, 15);
			
			extractedVerticesIDExtended.addAll(best_extension_path.subList(loc_in_list+1, best_extension_path.size()));
			//System.out.println("Path after " + extractedVerticesIDs.toString());
			//System.out.println("-------------------------------------------------------");
			
		}
	
		
		return extractedVerticesIDExtended;
	}
	
	
	
	private static HashMap<List<Integer>, Pair<Integer>> cuffMinPaths(final DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash, 
			DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis)
	{
		
		
		HashMap<List<Integer>, Pair<Integer>> transcripts = new HashMap<List<Integer>,Pair<Integer>>();
		
		
		
		Set<PairPath> pairPaths = new HashSet<PairPath>();
		Map<PairPath, Integer> pairPathToReadSupport = new HashMap<PairPath, Integer>();
		
	
		debugMes("Beginning cuffMinPaths",10);

		
		// populate pairPathToReadSupport: PairPath => readSupport
		// and pairPaths hashset: the list of all PairPaths
		
		
		populate_pairpaths_and_readsupport(combinedReadHash, pairPaths, pairPathToReadSupport);
		
		ArrayList<PairPath> pairPathsList = new ArrayList<PairPath>(pairPaths);
		
		Comparator<PairPath> pairPathOrderComparer = new Comparator<PairPath>() { // sort by first node depth in graph
			public int compare(PairPath a, PairPath b) {

				if (a.equals(b)) {
					return(0);
				}

				
				Integer a_index = a.getFirstID();
				Integer b_index = b.getFirstID();

				int f1 = getSeqVertex(graph, a_index)._node_depth;   // why using FinishingTime instead of DiscoveryTime?
				int f2 = getSeqVertex(graph, b_index)._node_depth;
				if( f1 < f2 )
					return -1;
				else if( f1 > f2 )
					return 1;

				// same node depth.
				if (a_index < b_index)
					return -1;
				else if (a_index > b_index)
					return 1;
				
				// same first node ID
				
				
				
				// check last node
				Integer a_last_index = a.getLastID();
				Integer b_last_index = b.getLastID();

				int l1 = getSeqVertex(graph,a_last_index)._node_depth;
				int l2 = getSeqVertex(graph,b_last_index)._node_depth;
				if (l1 < l2) {
					return(-1);
				}
				else if (l1 > l2) {
					return(1);
				}
				
				// same last node depth too.
				// compare their node identifiers
				if (a_last_index < b_last_index)
					return(-1);
				else if (a_last_index > b_last_index)
					return(1);
				
				
				// default
				// not the same paths, but same start node and last node DFS, so just order based on hashcode
				return ( (a.hashCode() < b.hashCode()) ? 1 : -1);
					
			}
		
		};
		
		
		Collections.sort(pairPathsList, pairPathOrderComparer);
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			debugMes("SORTED PAIRPATHS IN ORDER:", 15);
			for (PairPath p : pairPathsList) {
				debugMes("\t" + p, 15);
			}
			
			
		}
		
		
		HashSet<Integer> vertices = extract_vertex_list_from_PairPaths(pairPathsList);
		
		PairPath[] pairPathArr = pairPathsList.toArray(new PairPath[pairPathsList.size()]);
		
		HashSet<List<Integer>> cuff_input_paths = new HashSet<List<Integer>>();
		
		HashMap<PairPath,Integer> pairPathToIntVal = new HashMap<PairPath,Integer>();
		debugMes("All PairPaths sorted by DFS", 10);
		for (int i = 0; i < pairPathArr.length; i++) {
			pairPathToIntVal.put(pairPathArr[i], i);
			debugMes("CuffFly Input PairPath: " + pairPathArr[i] + "  <index: " + i + ">", 10);
			cuff_input_paths.add(pairPathArr[i].getPath1());
		}
		
		
		//0. remove containments
		List<Integer> containments = getContainments(pairPathArr);
		for (int i = 0; i < containments.size(); i++) {
			
			pairPathsList.remove(pairPathArr[containments.get(i)]);
		}
		
		// refresh after removing containments.
		pairPathArr = pairPathsList.toArray(new PairPath[pairPathsList.size()]);

		
		//1. build pair path graph
		int[][] dag = getPairPathDAG(graph, dijkstraDis, pairPathArr);
		
		//print dag
		debugMes("DAG------Before Removing Uncertainties ---------------------",10);
		for(int i = 0; i < dag.length; i++)
		{
			String dag_text = "";
			for(int j = 0; j < dag.length; j++)
			{
				dag_text += dag[i][j] + " ";
			}
			debugMes(dag_text, 10);
		}
		
		
		//2.1 remove uncertain pair paths
		ArrayList<Integer> uncertain = getUncertain(dag, pairPathArr);
		 
		debugMes("REMOVING UNCERTAINTIES: " + uncertain, 10);
		
		for(int i = 0; i < uncertain.size(); i++)
		{
		
			pairPathsList.remove(pairPathArr[uncertain.get(i)]);
		}
		
		HashSet<Integer> vertices_after_removed_uncertainties = extract_vertex_list_from_PairPaths(pairPathsList);
		
		if (vertices_after_removed_uncertainties.size() < vertices.size()) {
			int missing_node_count = vertices.size() - vertices_after_removed_uncertainties.size();
			debugMes("WARNING, MISSING: " + missing_node_count + " of " + vertices.size() + " nodes after removing uncertainties", 10);
			for (Integer v : vertices) {
				if (! vertices_after_removed_uncertainties.contains(v)) {
					debugMes("WARNING, MISSING NODE: After removing uncertainties, missing node from graph: " + v, 10);
				}
			}
		}
		
		
		// refresh again, after now removing the uncertain entries.
		pairPathArr = pairPathsList.toArray(new PairPath[pairPathsList.size()]);
		
		
		//print pair paths
		debugMes("PAIR PATHS------After Uncertain Removed--------------",10);
		for(int i = 0; i < pairPathsList.size(); i++)
		{
			debugMes("PairPathAfterUncertainRemoved "+ i + " " + pairPathArr[i].toString() + "  <index: " + pairPathToIntVal.get(pairPathArr[i]) + ">", 10);
		}
		
		// regenerate the dag now that the uncertain entries are removed.
		dag = getPairPathDAG(graph, dijkstraDis, pairPathArr);
		
		//print dag
		debugMes("DAG----After Uncertain Removed -----------------------",10);
		for(int i = 0; i < dag.length; i++)
		{
			String dag_text = "";
			for(int j = 0; j < dag.length; j++)
			{
				dag_text += dag[i][j] + " ";
			}
			debugMes(dag_text, 10);
		}
		
		//2.2 check transitivity
		if(!checkTransitivity(dag, pairPathArr, pairPathToIntVal))
		{
			throw(new RuntimeException("Graph is NOT transitive!"));
			
		}
		
		// remove self-matches
		for (int i = 0; i < dag.length; i++) {
			dag[i][i] = 0;
		}
		
		
		//2.3 get matching
		BipartiteMatching bp = new BipartiteMatching(dag.length, dag.length, dag);
		bp.maxMatching();
		
		int[] rightMatching = bp.getRightMatching();
		int[] leftMatching = bp.getLeftMatching();
		
		debugMes("Matching---------------------------------",10);
		bp.printRightMatching();
		
		
		//2.4 get chains from matching
		ArrayList<ArrayList<PairPath>> chains = new ArrayList<ArrayList<PairPath>>();
		boolean[] addedToChain = new boolean[rightMatching.length];
		Arrays.fill(addedToChain, false);
		ArrayList<PairPath> curChain;
		
		for (int i = 0; i < rightMatching.length; i++) {
			if (rightMatching[i] == -1) {
				// start new chain
				ArrayList<PairPath> chain = new ArrayList<PairPath>();
				chains.add(chain);
				chain.add(pairPathArr[i]);
				int j = i;
				while (leftMatching[j] != -1) {
					j = leftMatching[j];
					chain.add(pairPathArr[j]);
				}
			}
		}
		
		
		// report the chain info.
		HashMap<Integer,Boolean> seen = new HashMap<Integer,Boolean>();
		debugMes("Number of chains: " + chains.size(),10);
		for(int j = 0; j < chains.size(); j++)
		{
			curChain = chains.get(j);
			System.out.println("Chain: " + j);
			for(int k = 0; k < chains.get(j).size(); k++)
			{
				PairPath p = chains.get(j).get(k);
				Integer p_pos = pairPathToIntVal.get(p);
				String seen_text = (seen.containsKey(p_pos)) ? " *** ERROR, ALREADY INCLUDED IN ANOTHER CHAIN *** " : "";
				debugMes(p + " Pos:[" + p_pos + "] " + seen_text,10);
				seen.put(p_pos, new Boolean(true));
			}
			debugMes("",10);
		}
		
		// 3. foreach chain:
		// 3.1. extract nodes from chains
		DijkstraShortestPath dsp = new DijkstraShortestPath(graph);
		
		for(int i = 0; i < chains.size(); i++)
		{
			HashSet<Integer> extracted = new HashSet<Integer>();
			ArrayList<SeqVertex> extractedVertices = new ArrayList<SeqVertex>();
			
			List<List<Integer>> chain_i_path_list = new ArrayList<List<Integer>>();
			
			curChain = chains.get(i);
			for(int j = 0; j < curChain.size(); j++)
			{
				chain_i_path_list.add(curChain.get(j).getPath1());
				
				/*
				extracted.addAll(curChain.get(j).getPath1());
				if(curChain.get(j).hasSecondPath())
				    extracted.addAll(curChain.get(j).getPath2());
				    
				*/
					
			}
			
			List<List<Integer>> chain_i_collapsed_paths = Path.collapse_compatible_paths_to_min_set(chain_i_path_list);
			
			for (List<Integer> collapsed_path : chain_i_collapsed_paths) {

				transcripts.put(collapsed_path, new Pair(new Integer(1), new Integer(1)));
			
			}
			
			
			
			/*
			for(Integer id : extracted)
			{
				extractedVertices.add(getSeqVertex(graph, id));
			}
			//extractedVerticesIDs.addAll(extracted);
			
			// 3.2. sort according to topological order of BTFL graph 
			Collections.sort(extractedVertices, new SeqVertexFinishTimeComparator());
			
			String node_id_list_text = "";
			for (SeqVertex v : extractedVertices) {
				node_id_list_text += v.getID() + " ";
			}
			debugMes("Extracted vertices for chain: " + i + " and sorted is: " + node_id_list_text + "\n", 10);
			
			// Fill in any gaps
			// 3.3. path=[L_1]; For each i in 1:length(node_list)
			int j = 0;
			int num_vertices = extractedVertices.size(); // note, extractedVertices grows in size during iterations below.
			while(j < num_vertices - 1)
			{		
				
				//System.out.println(j);
				SeqVertex current = extractedVertices.get(j);
				SeqVertex next = extractedVertices.get(j + 1);
				
				// 3.3.1 if L_i == L_(i+1) then nothing
				// -There are no duplicates since extractedVertices was created from
				//  building the set of extracted vertex IDs
				
				// 3.3.2 else if exists an edge from L_i to L_(i+1) then nothing(?)
				if(graph.getSuccessors(current).contains(next)) {
					j++;
					continue;
				}
				
				// 3.3.3 else find a single path (p = L_i,..., L_(i+1)): 
				// append all P_j (j=2:end) to our path
				
				//List<SimpleEdge> sp = org.jgrapht.alg.DijkstraShortestPath.findPathBetween((Graph)graph, current, next);
				
				List<SimpleEdge> sp = dsp.getPath(current, next);
				debugMes("Found shorteset path between " + current.getID() + " and " + next.getID() + ":", 10);
				ArrayList<SeqVertex> toAdd = new ArrayList<SeqVertex>();
				for(SimpleEdge edge : sp) {
					SeqVertex v = graph.getDest(edge);
					toAdd.add(v);
					debugMes("\t" + v.getID(), 10);
				}
				toAdd.remove(next);
				extractedVertices.addAll(toAdd);
				
				j++;
			}
			
			ArrayList<Integer> extractedVerticesIDs = new ArrayList<Integer>();
			Collections.sort(extractedVertices, new SeqVertexFinishTimeComparator());
			for(SeqVertex v: extractedVertices)
			{
				//System.out.println("Adding vertex with ID: " + v.getID());
				extractedVerticesIDs.add(v.getID());
			}
			
			boolean extend_paths = ! CUFF_NO_EXTEND; // just for debugging purposes
			if (extend_paths) {
				ArrayList<Integer> extended = extendChain(extractedVerticesIDs, graph, pairPathToReadSupport, dijkstraDis);
				transcripts.put(extended, new Pair(new Integer(1), new Integer(1)));
			}
			else {
				transcripts.put(extractedVerticesIDs, new Pair(new Integer(1), new Integer(1)));
			}
			*/
			
			
			
		}	// end of foreach chain	
		
		
		
		debugMes("Cuff-based reconstructions of transcripts:", 10);
		for (List<Integer> p : transcripts.keySet()) {
			
			debugMes("CuffFly Output Path: " + p, 10);
			if (cuff_input_paths.contains(p))
				debugMes("\t** Original cuffpath",10);
			else
				debugMes("\t** NOT Original cuffpath", 10);
		}
		
		
		return transcripts;
	}
	
	
	
	private static HashSet<Integer> extract_vertex_list_from_PairPaths(
			ArrayList<PairPath> pairPathsList) {
		
		HashSet<Integer> vertices = new HashSet<Integer>();
		
		for (PairPath pp : pairPathsList) {
			


			for (List<Integer> path: pp.get_paths()) {

				for (Integer node : path) {
					vertices.add(node);
				}
			}
			
			
		}
		
		return(vertices);
		
		
	}


	private static List<Integer> getContainments(PairPath[] pairPathArr) {
		
		HashSet<Integer> containments = new HashSet<Integer>();
		
		for (int i = 0; i < pairPathArr.length; i++) {
			
			for (int j = 0; j < pairPathArr.length; j++) {
				
				if (i == j)
					continue;
				
				if (pairPathArr[i].isCompatibleAndContainedByPairPath(pairPathArr[j])) {
					debugMes("CONTAINMENT: " + pairPathArr[i] + " is contained by " + pairPathArr[j], 10);
					containments.add(i);
				}
				
				
			}
			
		}
		
		List<Integer> containment_list = new ArrayList(containments);
		
		return(containment_list);
		
		
	}


	private static void populate_pairpaths_and_readsupport(
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash,
			Set<PairPath> pairPaths,
			Map<PairPath, Integer> pairPathToReadSupport) {
		
		
		for(Integer i : combinedReadHash.keySet())
		{
			Set<PairPath> temp = combinedReadHash.get(i).keySet();
			for(PairPath j: temp)
			{
				if(j != null)
				{
					
					Integer read_count_support = combinedReadHash.get(i).get(j);
					
					j = j.trimSinkNodes();
					
					if (MAKE_PE_SE) { // separate the paths:
						
						// convert split paths into individual paths to avoid uncertainties.
						
						
						PairPath first_part = new PairPath(j.getPath1());
						if (! pairPaths.contains(first_part)) {
							pairPaths.add(first_part);
							pairPathToReadSupport.put(first_part, read_count_support);
						}
						else {
							// already there, just increment the read count support
							pairPathToReadSupport.put(first_part, pairPathToReadSupport.get(first_part)+read_count_support);
						}
						
						
						pairPathToReadSupport.put(first_part, read_count_support);
						
						if (j.hasSecondPath()) {
							PairPath second_part = new PairPath(j.getPath2());
							if (! pairPaths.contains(second_part)) {
								pairPaths.add(second_part);
								pairPathToReadSupport.put(second_part, read_count_support);
							}
							else {
								// already there, just increment the read count support
								pairPathToReadSupport.put(second_part, pairPathToReadSupport.get(second_part)+read_count_support);
							}
						}
						
						
						
					} else {
					
						// using pair paths instead of the split pairs (original)
						pairPaths.add(j);
						pairPathToReadSupport.put(j, read_count_support);
					
					}
				}
			}
		}
		
		
	}



	private static HashMap<List<Integer>, Pair<Integer>> butterfly (DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
		Set<SeqVertex> comp,
			HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash, 
			long totalNumReads,
			PrintStream pout_all,
			DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis, 
			DijkstraDistanceWoVer<SeqVertex,SimpleEdge> dijkstraDisWoVer,
			HashMap<Integer, List<List<Integer>>> tripletMapper,
			HashMap<Integer, List<List<Integer>>> extendedTripletMapper,
			HashMap<Integer,Boolean> xStructuresResolvedByTriplets
			
			)
	{



			
			Pair<HashMap<List<Integer>,Pair<Integer>>> FinalPathsPair = getAllProbablePaths(graph,comp,
					combinedReadHash,dijkstraDis,dijkstraDisWoVer, 
					tripletMapper,extendedTripletMapper, xStructuresResolvedByTriplets);

			HashMap<List<Integer>,Pair<Integer>> FinalPaths_diff = FinalPathsPair.getFirst();
			HashMap<List<Integer>,Pair<Integer>> FinalPaths_all = FinalPathsPair.getSecond();
			
			return FinalPaths_all;
	}



	/**
	 * given the graph, find all single nt bubbles, and choose the majority vote. 
	 * add the weights to the majority path, and add the prevID  
	 * v -> v1 -> vend
	 * v -> v2 -> vend
	 * @param graph
	 */
	private static void removeSingleNtBubbles(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {

		SeqVertex v1=null ,v2 = null, vend = null;
		SeqVertex vToKeep=null ,vToRemove = null;
		SimpleEdge e1ToKeep = null, e1ToRemove = null;
		SimpleEdge e2ToKeep = null, e2ToRemove = null;

		Vector<SeqVertex> removeV = new Vector<SeqVertex>();
		Collection<SeqVertex> allV = new HashSet<SeqVertex>();
		allV.addAll(graph.getVertices());

		for (SeqVertex v : allV)
		{
			if (removeV.contains(v))
				continue;
			if (graph.getSuccessorCount(v)==2)
			{
				Collection<SeqVertex> children = graph.getSuccessors(v);
				Iterator<SeqVertex> iter = children.iterator();
				v1 = iter.next();
				v2 = iter.next();


				int len1 = v1.getNameKmerAdj().length();
				int len2 = v2.getNameKmerAdj().length();

				debugMes("SNP_collapse candidates: " + v1  + " len: " + len1 + " and " + v2 + " len: " + len2, 15);
				
				
				if (len1==KMER_SIZE && len2==KMER_SIZE && 
						graph.getSuccessorCount(v1)==1 && 
						graph.getSuccessorCount(v2)==1 &&
						getSingleSuccessor(graph,v2).equals(getSingleSuccessor(graph,v1)))
				{
					vend = getSingleSuccessor(graph,v1);
					if (graph.findEdge(v, v1).getWeight() > graph.findEdge(v, v2).getWeight())
					{ //keep v1, loose v2
						vToKeep = v1;
						vToRemove = v2;
					}else
					{ //keep v2, loose v1
						vToKeep = v2;
						vToRemove = v1;
					}
					e1ToKeep = graph.findEdge(v, vToKeep);
					e2ToKeep = graph.findEdge(vToKeep, vend);
					e1ToRemove = graph.findEdge(v, vToRemove);
					e2ToRemove = graph.findEdge(vToRemove, vend);
					debugMes("SNP_collapse: merging the node "+vToRemove.getID()+" to the node "+vToKeep.getID(),15);

					SeqVertex newV = new SeqVertex(getNextID(), vToKeep.getName());
					newV.copyTheRest(vToKeep);
					newV.addToPrevIDs(vToKeep,vToRemove,LAST_REAL_ID);

					graph.addVertex(newV);

					graph.addEdge(new SimpleEdge(e1ToKeep.getWeight() + e1ToRemove.getWeight(), v.getID(), newV.getID()), v, newV);
					graph.addEdge(new SimpleEdge(e2ToKeep.getWeight() + e2ToRemove.getWeight(), newV.getID(), vend.getID()), newV,vend);

					removeV.add(vToRemove);
					removeV.add(vToKeep);


				}
			}
		}

		for (SeqVertex rv : removeV)
		{
			debugMes("removing the single nt variation vertex "+rv.getID(),20);
			graph.removeVertex(rv);
		}

	}




	/**
	 * given the graph, find all single nt bubbles, and choose the majority vote. 
	 * add the weights to the majority path, and add the prevID  
	 * v -> v1 -> vend
	 * v -> v2 -> vend
	 * @param graph
	 * @throws Exception 
	 */
	private static void removeSingleNtBubblesWithDegenerateCode(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) throws Exception {

		SeqVertex v1=null ,v2 = null, vend = null;
		SimpleEdge eTop1 = null, eTop2 = null;
		SimpleEdge eBottom1 = null, eBottom2 = null;

		Vector<SeqVertex> removeV = new Vector<SeqVertex>();
		Collection<SeqVertex> allV = new HashSet<SeqVertex>();
		allV.addAll(graph.getVertices());

		for (SeqVertex v : allV)
		{
			if (removeV.contains(v))
				continue;
			if (graph.getSuccessorCount(v)==2)
			{
				Collection<SeqVertex> children = graph.getSuccessors(v);
				Iterator<SeqVertex> iter = children.iterator();
				v1 = iter.next();
				v2 = iter.next();


				int len1 = v1.getName().length();
				int len2 = v2.getName().length();

				if (len1==1 && len2==1 && 
						graph.getSuccessorCount(v1)==1 && 
						graph.getSuccessorCount(v2)==1 &&
						getSingleSuccessor(graph,v2).equals(getSingleSuccessor(graph,v1)))
				{
					vend = getSingleSuccessor(graph,v1);
					String key;
					if (String.CASE_INSENSITIVE_ORDER.compare(v1.getName(),v2.getName())<0)
						key =  v1.getName()+v2.getName();
					else
						key =  v2.getName()+v1.getName();
					String name = getDegenerateRepresentation(key);
					SeqVertex newV = new SeqVertex(getNextID(), name);

					if (graph.findEdge(v, v1).getWeight() > graph.findEdge(v, v2).getWeight())
						newV.copyTheRest(v1);
					else
						newV.copyTheRest(v2);

					eTop1 = graph.findEdge(v, v1);
					eBottom1 = graph.findEdge(v1, vend);
					eTop2 = graph.findEdge(v, v2);
					eBottom2 = graph.findEdge(v2, vend);
					debugMes("merging the nodes "+v1.getID()+" and the node "+v2.getID()+" to the node "+newV,18);

					newV.addToPrevIDs(v1,v2,LAST_REAL_ID);
					newV.setFrequencies(v1.getName(),eTop1.getWeight(),v2.getName(),eTop2.getWeight());
					graph.addVertex(newV);
					graph.addEdge(new SimpleEdge(eTop1.getWeight() + eTop2.getWeight(), v.getID(), newV.getID()), v, newV);
					graph.addEdge(new SimpleEdge(eBottom1.getWeight() + eBottom2.getWeight(), newV.getID(), vend.getID()), newV,vend);

					removeV.add(v1);
					removeV.add(v2);


				}
			}
		}

		for (SeqVertex rv : removeV)
		{
			debugMes("removing the single nt variation vertex "+rv.getID(),20);
			graph.removeVertex(rv);
		}

	}

	/**
	 * return the single successor of this node in this graph
	 * @param graph
	 * @param v2
	 * @return
	 */
	private static SeqVertex getSingleSuccessor(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v) {
		Collection<SeqVertex> children = graph.getSuccessors(v);
		if (children.size()!=1)
			return null;

		SeqVertex vout = children.iterator().next();
		return vout;

	}


	/**
	 * find edges that are extremely high compared to both side (a single very abundant kmer, and fix their support
	 * @param graph
	 * @param inFlow 
	 * @param outFlow 
	 */
	private static void fixExtremelyHighSingleEdges(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, HashMap<Integer,Integer> outFlow, HashMap<Integer,Integer> inFlow) {

		debugMes("fixExtremelyHighSingleEdges()", 5);
		for (SimpleEdge e : graph.getEdges())
		{
			double supp =e.getWeight(); 
			Integer sourceID = graph.getSource(e).getID();
			Integer targetID = graph.getDest(e).getID();
			Integer inFlowToSource = inFlow.get(sourceID);
			Integer outFlowOfTarget = outFlow.get(targetID);

			if (inFlowToSource!= null && outFlowOfTarget!= null && 
					supp > inFlowToSource*EXTREME_EDGE_FLOW_FACTOR && supp > outFlowOfTarget*EXTREME_EDGE_FLOW_FACTOR)
			{
				double newSupp = Math.max(inFlowToSource, outFlowOfTarget);
				debugMes("the support of edge "+sourceID+"->"+targetID+" has changed from "+supp+" to "+newSupp,20);
				e.setWeight(newSupp);
			}

		}
	}


	/**
	 * given the graph and the final paths, find x structures that belong to only two paths, which resolve this structure.
	 * @param graph
	 * @param comp 
	 * @param finalPaths
	 * @return
	 */
	private static int countNumOfXstructuresResolved(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Set<SeqVertex> comp, HashMap<List<Integer>,Pair<Integer>> finalPaths) {

		int res = 0;

		for (SeqVertex v : comp)
		{
			if (graph.inDegree(v)>1 && graph.outDegree(v)>1)
			{
				//this is an x-structure
				int maxPaths = Math.max(graph.inDegree(v), graph.outDegree(v));
				
				Integer bef,after;
				int vid = v.getID();
				HashMap<Pair<Integer>,Integer> befAndAfterNodes = new HashMap<Pair<Integer>, Integer>();
				Pair<Integer> key;
				for (List<Integer> path : finalPaths.keySet())
				{
					int index = path.indexOf(vid);

					if (index!=-1 && index!=0 && index!=path.size()-1) // vid is not the first or the last
					{

						bef = path.get(index-1);
						after = path.get(index+1);
						key = new Pair<Integer>(bef,after);
						if (!befAndAfterNodes.containsKey(key))
							befAndAfterNodes.put(key,1);
						else
							befAndAfterNodes.put(key,befAndAfterNodes.get(key)+1);
					}
				}

				String triplets = "";
				for (Pair<Integer> befAndAfterNode : befAndAfterNodes.keySet()) {
					Integer before1 = (Integer) befAndAfterNode.getFirst();
					Integer after1 = (Integer) befAndAfterNode.getSecond();
					triplets += "[" + before1 + "-" + vid + "-" + after1 + "=" + befAndAfterNodes.get(befAndAfterNode) + "] "; 
				}
				
				if (befAndAfterNodes.keySet().size()==maxPaths)
				{
					debugMes("vertex "+v.getID()+" IS resolved in an X-structure: " + triplets, 10);
					res++;
				}
				else {
					debugMes("vertex " + v.getID() + " is NOT resolved in an X-structure: " + triplets, 10);
				}
			}
		}

		return res;
	}


	private static HashMap<Integer,Boolean> getXstructuresResolvedByTriplets (
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Set<SeqVertex> comp,
			HashMap<Integer, List<List<Integer>>> tripletMapper) {

		HashMap<Integer,Boolean> xStructuresResolvedByTriplets = new HashMap<Integer,Boolean>();

		for (SeqVertex v : comp)
		{
			if (graph.inDegree(v)>1 && graph.outDegree(v)>1)
			{
				//this is an x-structure
				Integer vertex_id = v.getID();
				
				if (tripletMapper.containsKey(vertex_id)) {
				
					debugMes("vertex " + vertex_id + " IS resolved in an X-structure: " + tripletMapper.get(vertex_id), 10);
					xStructuresResolvedByTriplets.put(vertex_id, true);
				}
				else {
					debugMes("vertex " + v.getID() + " is UN-resolved X-structure. ", 10);
					xStructuresResolvedByTriplets.put(vertex_id, false);
					
					
					
				}
			}
		}

		return (xStructuresResolvedByTriplets);
	}
	
	
	
	public static class FinalPaths implements Comparable<FinalPaths> {

		List<Integer> path;
		String sequence;

		public FinalPaths (List<Integer> p, String s) {
			path = p;
			sequence = s;
		}

		public int compareTo(FinalPaths f) {

			if (this.sequence.length() > f.sequence.length()) {
				return(-1);
			}
			else if (this.sequence.length() < f.sequence.length()) {
				return(1);
			}
			else {
				return(0);
			}

		}



	}

	

	
	

	/**
	 * Print all final paths
	 * @param finalPaths
	 * @param graph
	 * @param compID 
	 * @param p 
	 * @param name 
	 * @param totalNumReads 
	 * @param xStructuresResolvedByTriplets 
	 * @param separate_gene_ids 
	 * @throws FileNotFoundException 
	 */
	
	private static void printFinalPaths(
			HashMap<List<Integer>,Pair<Integer>> finalPaths,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			PrintStream p, 
			String name, 
			long totalNumReads,
			HashMap<List<Integer>,ArrayList<String>> final_paths_to_long_read_content, 
			HashMap<List<Integer>, Integer> separate_gene_ids) 

			throws FileNotFoundException 
		{

		

		debugMes("Final Paths: " + finalPaths.size(), 10);

		DecimalFormat df = new DecimalFormat("#.#");
		
	
		name = name.replace(".graph", "");

		
		HashMap<Integer,Integer> local_gene_id_mapping = new HashMap<Integer,Integer>();
		HashMap<Integer,Integer> local_seq_counter = new HashMap<Integer,Integer>();
		
		
		
		for (List<Integer> path : finalPaths.keySet()) {

			String seq = getPathSeq(graph,path);
			
			//print this path

			int gene_id;
			int iso_id;

			if (separate_gene_ids.containsKey(path)) {

				int local_gene_id = separate_gene_ids.get(path);

				if (local_gene_id_mapping.containsKey(local_gene_id)) {
					// seen this gene before
					gene_id = local_gene_id_mapping.get(local_gene_id);
					iso_id = local_seq_counter.get(gene_id) + 1;
					local_seq_counter.put(gene_id, iso_id);
				}
				else {
					gene_id = ++GENE_COUNTER;
					local_gene_id_mapping.put(local_gene_id, gene_id);
					iso_id = 1;
					// storing isoform id for the first time.
					local_seq_counter.put(gene_id, iso_id);
				}
			}
			else {
				// no gene clusters, each is unique:

				gene_id = ++GENE_COUNTER;
				iso_id = 1;
			}

			String seqName = name + "_g" + gene_id + "_i" + iso_id;

			String pathName = get_pathName_string(path, graph);

			seqName += " len="+seq.length() + " path="+ pathName;

			// Report the long read content information.
			if (final_paths_to_long_read_content.containsKey(path)) {
				// then got list of long read names assigned to this final path
				HashMap<PairPath, ArrayList<String>> long_read_paths_to_name_list = new HashMap<PairPath, ArrayList<String>>();


				for (String long_read_name : final_paths_to_long_read_content.get(path)) {

					PairPath pp = LONG_READ_NAME_TO_PPath.get(long_read_name);
					if (!long_read_paths_to_name_list.containsKey(pp)) {
						ArrayList<String> a = new ArrayList<String>();
						a.add(long_read_name);
						long_read_paths_to_name_list.put(pp, a);
					}
					else {
						long_read_paths_to_name_list.get(pp).add(long_read_name);
					}

				}
				seqName = seqName + " long_read_mappings: " + long_read_paths_to_name_list;
			}




			debugMes("\nFinal path reported: " + seqName, 10);
			p.print(getSeqFasta(seq, seqName));

		}


	}


	/*
	private static HashMap<List<Integer>, Boolean> remove_transcripts_with_insufficent_read_support(
			Set<List<Integer>> pathSet, HashMap<List<Integer>,Integer> seqLengthMap,
			HashMap<List<Integer>, Boolean> remove_low_expr_isoforms) {
	
		for (List<Integer> path : pathSet) {
			
			float path_frag_count = pc.get_transcript_to_sum_frag_counts(path);
			debugMes("PATH_TO_FRAG_COUNT: " + path_frag_count + ", FRAGS_PER_TRANS_LEN: " + path_frag_count/seqLengthMap.get(path)*100, 10);
			
		}
		
		
		return remove_low_expr_isoforms;
	}

	*/


	private static Object seqLengthMap() {
		// TODO Auto-generated method stub
		return null;
	}


	
	/*
	private static HashMap<List<Integer>, Float> get_pct_expr_isoforms_of_genes(
			HashMap<List<Integer>, Integer> separate_gene_ids) {
		
		// This relies on having run the EM to estimate relative expression.
		
		if (pc == null) {
			throw new RuntimeException("removal of low isoforms requires EM was run, but pc is null");
		}
		
		HashMap<List<Integer>, Float> pct_expr_isoform = new HashMap<List<Integer>,Float>();
		
		HashMap<Integer,Float> max_gene_expr_per_gene = new HashMap<Integer,Float>();
		
		for (List<Integer> transcript : separate_gene_ids.keySet()) {
			Integer gene_id = separate_gene_ids.get(transcript);
			
			Float expr = pc.get_expr(transcript);
			if (max_gene_expr_per_gene.containsKey(gene_id)) {
				if (max_gene_expr_per_gene.get(gene_id) < expr)
					max_gene_expr_per_gene.put(gene_id, expr);
			}
			else
				max_gene_expr_per_gene.put(gene_id, expr);
			
		}
		
		for (List<Integer> transcript : separate_gene_ids.keySet()) {
			Integer gene_id = separate_gene_ids.get(transcript);
			
			Float expr = pc.get_expr(transcript);
			Float max_gene_expr = max_gene_expr_per_gene.get(gene_id);
			
			float pct_isoform_expr = expr/max_gene_expr * 100;
			
			debugMes("Relative expression: " + pct_isoform_expr + ", gene: " + gene_id + ", path: " + transcript, 10);
		
			
			pct_expr_isoform.put(transcript, pct_isoform_expr);
			
		}
		return(pct_expr_isoform);
		
	}
	*/



	private static HashMap<List<Integer>, Integer> group_paths_into_genes(
			HashMap<List<Integer>, Pair<Integer>> finalPaths_all, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
	
		debugMes("Grouping paths into genes", 10);

		HashMap<Integer,Integer> node_length_map = new HashMap<Integer,Integer>(); // track node lengths for seq pair comparisons
		HashMap<List<Integer>, Integer> seqLengthMap = new HashMap<List<Integer>,Integer>();
		Vector<FinalPaths> path_vec = new Vector<FinalPaths>();
		
		// get node lengths
		for (List<Integer> path : finalPaths_all.keySet())
		{
			String seq = getPathSeq(graph,path);

			seqLengthMap.put(path, seq.length());
			//System.out.println(seq);

			FinalPaths f = new FinalPaths(path, seq);
			path_vec.add(f);
			
			for (Integer nodeID : path) {
				if (nodeID > 0) {
					int node_length = getSeqVertex(graph, nodeID).getName().length();
					node_length_map.put(nodeID, node_length);
				}
				else 
					node_length_map.put(nodeID, 0);
			}
		}

		
		UndirectedSparseGraph<List<Integer>, String> sparseGraph = new UndirectedSparseGraph<List<Integer>, String>();
		
		List<List<Integer>> paths = new ArrayList<List<Integer>>(finalPaths_all.keySet());
		
		for (int i = 0; i <= paths.size()-2; i++) {
			
			List<Integer> path_i = paths.get(i);
			
			if (! sparseGraph.containsVertex(path_i))
				sparseGraph.addVertex(path_i);
				
			int path_i_len = 0;
			for (Integer node : path_i) {
				path_i_len += node_length_map.get(node);
			}
			for (int j = i + 1; j <= paths.size()-1; j++) {
				
				List<Integer> path_j = paths.get(j);
				if (! sparseGraph.containsVertex(path_j))
					sparseGraph.addVertex(path_j);
				
				
				int path_j_len = 0;
				
				int nodes_same_length = 0;
				for (Integer node : path_j) {
 					path_j_len += node_length_map.get(node);
				
 					if (path_i.contains(node))
 						nodes_same_length += node_length_map.get(node);
				}
	
				
				float iso_pct_overlap = Math.max((float)nodes_same_length / path_i_len * 100,
						(float)nodes_same_length / path_j_len * 100);
				
				debugMes("Isoform_overlap: Path_i:" + path_i + ", Path_j: " + path_j + ", overlap = " + iso_pct_overlap + "%", 10);
				
				
				if ( iso_pct_overlap >= MIN_ISOFORM_PCT_LEN_OVERLAP) {
					
					sparseGraph.addEdge("e_" + i + "_" + j, path_i, path_j);
					debugMes("IsoformEdge linking: " + path_i + " to " + path_j, 10);
					
				}
				
				
			}
		}
		
		HashMap<List<Integer>, Integer> gene_grouping = new HashMap<List<Integer>, Integer>();
		
		WeakComponentClusterer<List<Integer>, String> compClus = new WeakComponentClusterer<List<Integer>,String>();
		Set<Set<List<Integer>>> comps = compClus.transform(sparseGraph);
		
		debugMes("IsoformClustering, number of clusters = " + comps.size(), 10);
		
		// add the singletons back in
		HashMap<List<Integer>,Boolean> inCluster = new HashMap<List<Integer>,Boolean>();
				
		if (comps.size() == 0) {
			// all related
			for (List<Integer> path : finalPaths_all.keySet()) {
				gene_grouping.put(path, 1);
			}
			return(gene_grouping);
		}
		
		
		int cluster_count = 0;
		for (Set<List<Integer>> cluster : comps) {
			cluster_count++;
			for (List<Integer> path : cluster) {
				gene_grouping.put(path, cluster_count);
				debugMes("GeneCluster[" + cluster_count + "] contains: " + path, 10);
				inCluster.put(path, true);
			}
		}
		for (List<Integer> path : paths) {
			if (! inCluster.containsKey(path)) {
				cluster_count++;
				gene_grouping.put(path, cluster_count);
			}
		}
			
		
		return(gene_grouping);
		
	}
	

	/**
	 * given a path in the graph, return its sequence
	 * @param graph
	 * @param path
	 * @return
	 */
	private static String getPathSeq(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, List<Integer> path) {
		
		String seq = "";

		Boolean first_node = true;

		for (Integer nodeID : path) {
			if (nodeID>=0) {
				String node_seq = getSeqVertex(graph, nodeID).getName();
				if (! first_node) {
					node_seq = node_seq.substring(SeqVertex.get_kmer_length() -1);
						
				}
				first_node = false;
				seq += node_seq;
				
				//System.out.println("Node: " + nodeID + " has seq: " + node_seq);
			}	
			
		}
		

		return seq;
	}




	/**
	 * For each path of a read pair, ask how many reads support it.
	 * @param graph
	 * @param readNameHash
	 * @param dijkstraDis 
	 * @return
	 */
	private static HashMap<Integer, HashMap<PairPath, Integer>> getSuffStats_wPairs(
			DirectedSparseGraph<SeqVertex,SimpleEdge> graph, 
			HashMap<String, List<Read>> readNameHash, 
			DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis) {
		
		HashMap<Integer,HashMap<PairPath,Integer>>   combinedReadHash = new HashMap<Integer,HashMap<PairPath,Integer>>  ();

		Set<String> usedReads = new HashSet<String>();
		List<Read> curList = null;
		int numReadsUsed = 0;

		int numSingletons = 0;
		int numPairs = 0;
		int numPairsDiscarded = 0;

		for (String name : readNameHash.keySet())
		{
			if (usedReads.contains(name))
				continue; // ignoring reduncancy in the read set?

			curList = readNameHash.get(name);
			if (curList.size()==1)
			{//single read
				
				// *****************************
				// **  Single Read Processing **
				// *****************************
				
				Read r = curList.get(0);
				PairPath path = new PairPath(r.getPathIDs());
				Integer firstV = path.getFirstID();

				if (!combinedReadHash.containsKey(firstV))
					combinedReadHash.put(firstV, new HashMap<PairPath,Integer>()); // init

				if (!combinedReadHash.get(firstV).containsKey(path))
					combinedReadHash.get(firstV).put(path, 0);

				Integer counts = combinedReadHash.get(firstV).get(path);
				combinedReadHash.get(firstV).put(path,++counts); // count read having this path
				numReadsUsed++;
				debugMes("we have "+combinedReadHash.get(firstV).get(path)+" reads supporting the path: "+path,19);
				numSingletons++;

				// examine for long read.
				if (name.startsWith("LR$|")) { //   r.getSeq().length() >= MIN_LONG_READ_LENGTH) {
					LONG_READ_NAME_TO_PPath.put(name,  path);
					debugMes("LONG_READ_IDENTIFIED: " + name + " , path: " + path, 12);
					if (! LONG_READ_PATH_MAP.containsKey(path)) {
						ArrayList<String> nameList = new ArrayList<String>();
						nameList.add(name);
						LONG_READ_PATH_MAP.put(path, nameList);
					}
					else {  
						ArrayList<String> nameList = (ArrayList<String>) LONG_READ_PATH_MAP.get(path);
						nameList.add(name);
					}
				}


			}else { // paired read
				
				// *****************************
				// **  Paired Read Processing **
				// *****************************

				
				Read r1 = curList.get(0);
				List<Integer> path1 = r1.getPathIDs();

				Read r2 = curList.get(1);
				List<Integer> path2 = r2.getPathIDs();

				PairPath combinedPath = new PairPath(path1, path2);
				
				/*  move this to after repeat unrolling.
				 * for now, just keep it simple, store paths.
				 * 
				PairPath  combinedPath = combinePaths(graph,path1,path2,dijkstraDis);
				if (combinedPath.isEmpty())
				{
					debugMes("the paths "+path1+" and "+path2+" couldn't be combined",15);
					numPairsDiscarded++;
					continue;
				}
				*/


				Integer firstV = combinedPath.getFirstID();

				if (!combinedReadHash.containsKey(firstV))
					combinedReadHash.put(firstV, new HashMap<PairPath,Integer>()); //init

				if (!combinedReadHash.get(firstV).containsKey(combinedPath))
					combinedReadHash.get(firstV).put(combinedPath, 0); //add pairpath

				Integer counts = combinedReadHash.get(firstV).get(combinedPath);
				combinedReadHash.get(firstV).put(combinedPath,++counts); // increment counts for pairpath
				debugMes("we have "+combinedReadHash.get(firstV).get(combinedPath)+" reads supporting the path: "+combinedPath,18);

				numReadsUsed++;
				numPairs++;


			}
			usedReads.add(name);
		}
		debugMes("number of reads used = "+numReadsUsed,15);

		debugMes("## Read PathPair results: " + numSingletons + " singletons, "
				+ " num pairs: " + numPairs + ", num pairs discarded: " + numPairsDiscarded, 10);

		return combinedReadHash;
	}

	/**
	 * Given the graph, and two paths of the two reads, combine them into a single path
	 * @param graph
	 * @param path1
	 * @param path2
	 * @param dijkstraDis
	 * @return
	 */
	private static PairPath combinePaths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			List<Integer> path1, List<Integer> path2, DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis) {
		
		
		debugMes("combinePaths: " + path1 + ", " + path2, 15);
		
		
		SeqVertex firstV1 = getSeqVertex(graph,path1.get(0));
		SeqVertex lastV1 = getSeqVertex(graph,path1.get(path1.size()-1));
		SeqVertex firstV2 = getSeqVertex(graph,path2.get(0));
		SeqVertex lastV2 = getSeqVertex(graph,path2.get(path2.size()-1));
		PairPath path  = new PairPath();

		if (path1.containsAll(path2))
			path.setPath1(path1);
		else if (path2.containsAll(path1))
			path.setPath2(path2);  // note, gets moved to path1 later.
		
		//path1 --> path2
		else if (SeqVertex.isAncestral(lastV1, firstV2,dijkstraDis)>0 
				&&
				! lastV1.equals(firstV2)
			)
		{
			path.setPath1(path1);
			path.setPath2(path2);
		}
		//path2 --> path1
		else if (SeqVertex.isAncestral(lastV2, firstV1,dijkstraDis)>0
				&&
				! lastV2.equals(firstV1)
			)
		{
			path.setPath1(path2);
			path.setPath2(path1);
		}

		else if (SeqVertex.isAncestral(firstV2,firstV1,dijkstraDis)==0 && 
				SeqVertex.isAncestral(lastV2,lastV1,dijkstraDis)==0)
		{
			//there is no consistent path between read1 and read2
		}

		// ******************
		// Overlapping paths
		// ******************
		
		//path1(partial) -> path2
		else if (SeqVertex.isAncestral(firstV1,firstV2,dijkstraDis)>0 && 
				path1.indexOf(firstV2.getID())>=0)
		{
			int i = path1.indexOf(firstV2.getID());
			path.setPath1(path1.subList(0, i));
			path.addToPath1(path2);
		}

		//path2(partial) -> path1
		else if (SeqVertex.isAncestral(firstV2,firstV1,dijkstraDis)>0 &&
				path2.indexOf(firstV1.getID())>=0)
		{
			int i = path2.indexOf(firstV1.getID());
			path.setPath1(path2.subList(0, i));
			path.addToPath1(path1);
		}

		if (path.getPath1().isEmpty() && !path.getPath2().isEmpty())
			path.movePath2To1();

		
		// Try to impute connecting paths from path 1 to path 2
		if ((! path.getPath1().isEmpty()) && (! path.getPath2().isEmpty())) {
			SeqVertex fV1 = getSeqVertex(graph,path.getPath1().get(0));
			SeqVertex lV1 = getSeqVertex(graph,path.getPath1().get(path.getPath1().size()-1));
			SeqVertex fV2 = getSeqVertex(graph,path.getPath2().get(0));
			SeqVertex lV2 = getSeqVertex(graph,path.getPath2().get(path.getPath2().size()-1));
			
			debugMes("Examining imputation of path connecting pairpath: " + path + " nodes "+ lV1.getID() + " to " + fV2.getID(), 20);
			
			if (SeqVertex.isAncestral(lV1, fV2, dijkstraDis) > 0) {
				// note could return false if have a sequencing gap 
				
				boolean canExtend = true;
				SeqVertex v = lV1;
				// walk towards fV2
				List<Integer> intervening_vertex_ids = new ArrayList<Integer>();
				boolean impute = true;
				while (canExtend) {
					SeqVertex next = null;
					int count_connectable = 0;
					for (SeqVertex successor : graph.getSuccessors(v)) {
						if (SeqVertex.isAncestral(successor, fV2, dijkstraDis) > 0) {
							count_connectable++;
							next = successor;
						}
					}
					if (next != null && count_connectable == 1) {
						if (fV2.equals(next)) {
							// reached fV2
							break;
						}
						else {
							intervening_vertex_ids.add(next.getID());
						}
						v = next;

					}
					else {
						// either no connection or too many potential connections
						canExtend = false;
						impute = false;
					}

				}
				if (impute) {
					debugMes("Could Impute path connecting" + path + " containing intervening nodes: " + intervening_vertex_ids, 16);
					if (! intervening_vertex_ids.isEmpty()) {
						path.getPath1().addAll(intervening_vertex_ids);
					}
					path.getPath1().addAll(path.getPath2());
					path.getPath2().clear();

				}
				else {
					debugMes("Could not impute intervening nodes", 20);
				}
			}
			
		}
		
		
		return path;

	}



	/**
	 * Count how many vertices we have with in degree >1 & out degree >1
	 * @param graph
	 * @return
	 */
	private static int countNumOfXstructures(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		int res = 0;

		for (SeqVertex v : graph.getVertices())
		{
			if (graph.inDegree(v)>1 && graph.outDegree(v)>1)
				res++;
		}

		return res;
	}

	private static void printPairPaths(HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash){
		printPairPaths(combinedReadHash, "PAIRPATH");

	}

	private static void printPairPaths(HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash,
				String out_token){
		
		for(Map.Entry<Integer, HashMap<PairPath, Integer>> entry : combinedReadHash.entrySet())
		{
			debugMes("Start Vertex:" + entry.getKey(), 10);
			for(Map.Entry<PairPath, Integer> paths : entry.getValue().entrySet())
			{
				//System.out.println(entry.getValue());
				debugMes(out_token + ": " + paths, 10);
			}
		}
		
	}
	
	/**
	 * Given the graph and the hash with all reads, find all probable paths from S to T.
	 * @param graph
	 * @param comp 
	 * @param combinedReadHash
	 * @param dijkstraDis 
	 * @param dijkstraDisWoVer 
	 */
	private static Pair<HashMap<List<Integer>, Pair<Integer>>> getAllProbablePaths (
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Set<SeqVertex> comp, 
			HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash, 
			DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis, 
			DijkstraDistanceWoVer<SeqVertex,SimpleEdge> dijkstraDisWoVer,
			HashMap<Integer, List<List<Integer>>> tripletMapper,
			HashMap<Integer, List<List<Integer>>> extendedTripletMapper,
			HashMap<Integer,Boolean> xStructuresResolvedByTriplets
			) {

		debugMes("\nSECTION\n###############\n## Starting Butterfly Assembly ##\n###################\n\n", 5);
		
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			debugMes("PairPaths to assemble:", 15);
			printPairPaths(combinedReadHash, "PairPaths@BflyStart");
		}
		
		
		// paths that are constructed by tracing paths of reads through the graph
		HashMap<SeqVertex,List<List<Integer>>> Paths = new HashMap<SeqVertex,List<List<Integer>>>();
		
		// this holds reads of path until V + reads starting at V
		HashMap<List<Integer>,HashMap<PairPath,Integer>> PathReads = new HashMap<List<Integer>,HashMap<PairPath,Integer>>();
		
		HashMap<List<Integer>, HashSet<PairPath>> PathContainedReads = new HashMap<List<Integer>,HashSet<PairPath>>();
		
		
		HashMap<List<Integer>,Boolean> Extensions = new HashMap<List<Integer>,Boolean>();

		
		// final paths to be captured and reported.
		HashMap<List<Integer>,Pair<Integer>> FinalPaths_diff = new HashMap<List<Integer>,Pair<Integer>>();
		HashMap<List<Integer>,Pair<Integer>> FinalPaths_all = new HashMap<List<Integer>,Pair<Integer>>();

		/*
		ROOT.setDFS_FinishingTime(Integer.MAX_VALUE);
		T_VERTEX.setDFS_FinishingTime(-1);
	    */


		ROOT.setDepth(-1);
		T_VERTEX.setDepth(Integer.MAX_VALUE);
		
		/*
		SeqVertexFinishTimeComparator finishingTimeComparator = new SeqVertexFinishTimeComparator();

		PriorityQueue<SeqVertex> C = new PriorityQueue<SeqVertex>(comp.size(),finishingTimeComparator);

		*/
		
		PriorityQueue<SeqVertex> BflyQueue = new PriorityQueue<SeqVertex>(comp.size(), new SeqVertexNodeDepthComparator());

		
		// start the path listing out with a path containing the ROOT node only.
		BflyQueue.add(ROOT);
		List<Integer> tmpL = new ArrayList<Integer>();
		tmpL.add(ROOT.getID());
		ArrayList<List<Integer>> tmpPathList = new ArrayList<List<Integer>>();
		tmpPathList.add(tmpL);
		Paths.put(ROOT, tmpPathList);
		SeqVertex v;


		int total_num_nodes = comp.size();
		
		// *******************************
		// beginning path constructions
		// *******************************
		
		
		HashMap<Integer,Boolean> node_ID_visited = new HashMap<Integer,Boolean>();
		
		
		int num_nodes = 0;
		
		String Crep;
		while (!BflyQueue.isEmpty())
		{
	
			
			/*
			if (BFLY_GLOBALS.VERBOSE_LEVEL>=20)
			{
				Crep = "[";
				for (SeqVertex vp : C)
					Crep = Crep + "" +vp.getID()+":"+vp.getDFS_FinishingTime()+",";
				Crep += "]";
				debugMes("C = "+Crep,10);
			}
			*/
			
			debugMes("QUEUE IS: " + BflyQueue, 12);
			
			v = BflyQueue.poll(); 

			if (v.getID() > 0) {
				// if it has successors that haven't been visited yet, delay targeting it.
				List<SeqVertex> delay_tackle_vertices = new ArrayList<SeqVertex>();
				while ( (! BflyQueue.isEmpty()) &&  (! parents_all_visited(v, node_ID_visited, graph) ) ) {
					debugMes("* delaying tackling vertex: " + v.getID() + " since a parent hasn't been visited yet.", 12);
					delay_tackle_vertices.add(v);
					v = BflyQueue.poll();
				}
				if (BflyQueue.isEmpty() && ! parents_all_visited(v, node_ID_visited, graph)) {
					throw new RuntimeException("ERROR, queue ran out of nodes and current node has unvisited parents.");
				}
				if (! delay_tackle_vertices.isEmpty()) {
					// add them back to the queue
					BflyQueue.addAll(delay_tackle_vertices);
				}
			}

			
			// track the nodes we visit, avoid looping by extending from a node encountered earlier.  Loops should be handled long before here.
			if (node_ID_visited.containsKey(v.getID())) {
				debugMes("** already visited node in queue: " + v.getID(), 5);
				continue;
			}
			else {
				node_ID_visited.put(v.getID(), true);
			}

			debugMes("\n\n#### getAllProbablePaths() The next node in the queue C is "+v.getID(),5);
			
			num_nodes++;
			
			float pct_done = (float) num_nodes / total_num_nodes * 100;
			debugMes("\tbutterfly pct done: " + num_nodes + " / " + total_num_nodes + " = "  + pct_done + "% pct done.", 5);
			
			
			// get read paths that start at vertex V
			HashMap<PairPath,Integer> readsStartingAtV = combinedReadHash.get(v.getID());
			
			
			if (readsStartingAtV == null) {
				debugMes("ReadsStartingAtV_START_BFLY" + v.getID() +  " EMPTY", 15);
			}
			else {
				for (PairPath read : readsStartingAtV.keySet()) {
					debugMes("ReadsStartingAtV_START_BFLY, Node: " + v.getID() + " read: " + read, 15);
				}
			}
			
			/////////////////////////////////////////////////////////////
			// prep data structures required.
			// go over all paths of P[v], add all reads that start at v
			//////////////////////////////////////////////////////////////
			
			debugMes("Exploring extension of: " + Paths.get(v).size() + " paths that end at vertex: " + v.getID(), 5);
			
			
			//describe paths:
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
				debugMes("\n== Current Paths Constructed Up To Vertex: " + v.getID() + " :", 15);
				for (List<Integer> path: Paths.get(v)) {
					debugMes("PathPartialReconstruction@[" + v.getID() + "] : " + path, 15); 
					
				}
			}
			
			for (List<Integer> path : Paths.get(v))
			{
				if (!PathReads.containsKey(path))
					PathReads.put(path, new HashMap<PairPath,Integer>()); // init

				if (! PathContainedReads.containsKey(path)) 
					PathContainedReads.put(path, new HashSet<PairPath>());
				
				if (readsStartingAtV!=null && !readsStartingAtV.isEmpty())
				{
					debugMes("\nAdding the reads " +readsStartingAtV +" to the path "+ path, 17);
					PathReads.get(path).putAll(readsStartingAtV);
				
					/*
					// verify:
					for (PairPath pp : readsStartingAtV.keySet()) {
						debugMes("VERIFYING: " + v.getID() + " ReadStartingAtV: " + pp + " = " + PathReads.get(path).get(pp), 10);
					}
					*/
					
					// path that ends at V is associated with all reads that start at V
				}

				//keep track of all extensions
				Extensions.put(path, false);

				/*
				for (PairPath read : PathReads.get(path).keySet()) {
					debugMes("PATH: " + path + " initially stocked with read: " + read, 10);
				} 
				*/
				
				
			}


			//////////////////////////////////////////////////////
			//  Examine each path, try to extend by successor u
			/////////////////////////////////////////////////////
			

			// go over all descendants of v
			for (SeqVertex u : graph.getSuccessors(v))
			{

				debugMes("\n\n################################################\n"
						 + "###### Exploring extension of v: " + v.getID() + " by successor: " + u.getID() 
						 +"\n################################################\n\n", 12);
				
				if (! (comp.contains(u) || u.equals(T_VERTEX))) {
					debugMes("component either lacks: " + u.getID() + " or at sink", 12);
					continue; // only examine successor vertices that are contained within this subcomponent
				}
				
				int path_counter = 0;

				boolean vExtendedToU = false;

				List<List<Integer>> paths_ending_at_v = new ArrayList<List<Integer>>(Paths.get(v));
				debugMes("Count of paths ending at v: " + v.getID() +  " = " + paths_ending_at_v.size(), 12);
				
				// sort paths by pair-path support descendingly
				PathReadSupportComparator local_pc = new PathReadSupportComparator(PathReads);
				Collections.sort(paths_ending_at_v, local_pc);
				Collections.reverse(paths_ending_at_v); // now descending according to read support.
				
				
				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
					for (List<Integer> path : paths_ending_at_v) 
					{
						debugMes("path_ending_at_v: " + path, 15);
					}
					
				}
				
				
				// Examine each of the growing paths that end at vertex V for extension.
				// extend V by U if:
				//     a.  in TRIPLET LOCK mode and there exists a read-supported structure w-v-u
				//     b.  extended triplet exists and the growing path is consistent with the
				//         complex read path consistent with  [a.b.c]..w-v-u
				//     c.  haven't reached the limit of number of paths artificially capped  //FIXME: re-examine this.
				//     d.  ALL POSSIBLE PATHS parameter set
				//     e.  path has enough read support: look back the number of nodes before U in this growing path
				//         that yield at least the required path support sequence distance, and ensure read paths
				//         demonstrate compatibility/containment
				
				for (List<Integer> path : paths_ending_at_v) 
				{

					// remember, only looking at extensions that contain 'u' now.

					debugMes("\n\n# [PathCounter(" + u.getID() + ")=" + path_counter + " Examining potential extension of path ending at node V: " + v.getID() 
							+ " by successor: " + u.getID() 
							+ ", via path=" + path, 15);
					
					Boolean path_wvu_acceptable = true; // by default
					Boolean extended_triplet_path_compatible = false;

					if (path.size() >= 3) {
						Integer w = path.get(path.size()-2); // create triplet w-v-u  

						if (tripletMapper.containsKey(v.getID()) 
								&& tripletMapper.get(v.getID()).size() > 1)  // at least partially resolved structure via read path
						{
							
							List<Integer> triplet = new ArrayList();
							triplet.add(w); // left
							triplet.add(v.getID()); // central
							triplet.add(u.getID()); // right
							List<List<Integer>> triplet_list = tripletMapper.get(v.getID());
							if (tripletSupported(triplet_list, triplet)){
								// Hurray, got triplet support
								debugMes("Triplet Path: " + triplet + " *IS* supported by reads.", 15);
								path_wvu_acceptable = true;
								
								// do extended triplet search.
								// path must be compatible with at least one of the complex prefix paths
								// ensuring compatible with the larger path context that may extend beyond a triplet
								
								List<Integer> pathWu = new ArrayList<Integer>();  // pathWu = path with u
								pathWu.addAll(path);
								pathWu.add(u.getID());

								//////////////////////////
								// extended triplet search
								//////////////////////////
								
								for (List<Integer> prefix_path : extendedTripletMapper.get(u.getID())) {
									PairPath ppath = new PairPath(prefix_path);

									if (ppath.isCompatibleAndContainedBySinglePath(pathWu)) {

										debugMes("EXTENDED_TRIPLET_SEARCH: " + ppath + " compared to " + pathWu + " True", 15);
										extended_triplet_path_compatible = true;
										break;
									}
									else {
										debugMes("EXTENDED_TRIPLET_SEARCH: " + ppath + " compared to " + pathWu + " False", 15);	
									}

								}
								
							}
							else {
								// lock down node, don't allow alternative structures not supported by reads here.
								debugMes("Triplet Path: " + triplet + " is *NOT* supported by reads.", 15);
								path_wvu_acceptable = false;
							}
						}
						else {
							debugMes("TripletMapper doesnt contain node: " + v.getID(), 15);
							// if node v is at center of X-structure and there are no valid triplets, disable extension
							if (FRACTURE_UNRESOLVED_XSTRUCTURE && xStructuresResolvedByTriplets.containsKey(v.getID())) {
								debugMes("Node " + v.getID() + " is at center of X structure and no triplet support detected.  FractureUnresolvedX set, so Disabling extension.", 10);
								path_wvu_acceptable = false;
							}
						}
					}
					else {
						debugMes("path " + path + " is too short to check for triplet support.", 15);
					}


					HashMap<PairPath,Integer> readsOfPathUntilV = PathReads.get(path); //this holds reads of path until V + reads starting at V

					if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
						debugMes("ReadsOfPathUntilV: PATH: " + path, 15);
						for (PairPath pp : readsOfPathUntilV.keySet())
							debugMes("ReadsOfPathUntiV: READ: " + pp, 15);
					}
					
					if (  
							path_wvu_acceptable // path_wvu only matters under triplet locking mode.

							&&

							// place restriction on the number of paths pursued here.
							// if extended triplet compatible, be sure to explore it... dont want to lose a good path here.
							(extended_triplet_path_compatible || path_counter <= MAX_NUM_PATHS_PER_NODE_EXTEND)  

							&&

							(
									ALL_POSSIBLE_PATHS 
									|| 	
									pathHasEnoughReadSupport(readsOfPathUntilV,path,u,graph,dijkstraDisWoVer)
									||
									u.getID() < 0 // a sink node, if path made it this far, sink can be added.
									)



							)

					{

						path_counter++;

						// add [path,u] to paths of u.  Each vertex contains all the paths that led up to it.
						if (!Paths.containsKey(u))
							Paths.put(u, new ArrayList<List<Integer>>());

						List<Integer> pathWu = new ArrayList<Integer>();  // pathWu = path with u
						pathWu.addAll(path);
						pathWu.add(u.getID());
						if (!Paths.get(u).contains(pathWu)){
							debugMes("\nSuccessful extension of " + u.getID() + " to generate path " +pathWu, 15);
							Paths.get(u).add(pathWu);
						}

						//update reads of [path,u] : includes all reads that are consistent with path Wu up to and including u
						updateReadsOfPath(PathReads,PathContainedReads, pathWu,readsOfPathUntilV,u.getID(),graph,dijkstraDis);

						//update extension
						Extensions.put(path, true);
						vExtendedToU = true;
					}
					else {
						debugMes("No extension of path " + path + " by " + u, 15);
						if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
							boolean pathHasEnoughSupport = pathHasEnoughReadSupport(readsOfPathUntilV,path,u,graph,dijkstraDisWoVer);
									debugMes("\tpath_counter = " + path_counter + ", path_wvu_acceptable=" + path_wvu_acceptable
											+ " pathHasEnoughReadSupport=" + pathHasEnoughSupport, 15);
						}
						
					}

				}  // end of extending paths that end at V by U
				
				
				
				if (!BflyQueue.contains(u))
				{
					debugMes(u.getID()+" was added to the queue",17);
					BflyQueue.add(u);

				}
				//if v didn't extend to u, and we have an edge there, add (v,u) as a new path
				if ( (!vExtendedToU) )
				{
					debugMes("the edge (v-u) was not used in any extension: "+v.getID()+"->"+u.getID(),15);
					if (!Paths.containsKey(u))
						Paths.put(u, new ArrayList<List<Integer>>());
					List<Integer> vuPath = new ArrayList<Integer>();
					vuPath.add(v.getID());
					vuPath.add(u.getID());

					Paths.get(u).add(vuPath);

					//add the reads
					if (!PathReads.containsKey(vuPath))
						PathReads.put(vuPath, new HashMap<PairPath,Integer>());

					if (readsStartingAtV!=null && !readsStartingAtV.isEmpty())
					{
						debugMes("adding the reads " +readsStartingAtV +" to the path "+ vuPath, 17);
						PathReads.get(vuPath).putAll(readsStartingAtV);
						updateReadsOfPath(PathReads, PathContainedReads, vuPath,readsStartingAtV,u.getID(),graph,dijkstraDis);

					}

				}

			}  // end of exploration of successors U of V
			
			
			
			//report the paths that were not extended AND remove them from Paths
			List<List<Integer>> removePaths = new ArrayList<List<Integer>>();
			for (List<Integer> path : Paths.get(v))
			{
				SeqVertex lastV = getSeqVertex(graph, path.get(path.size()-1));

				if (!lastV.equals(T_VERTEX) && Extensions.get(path)!=null && !Extensions.get(path))
				{
					if (getSeqPathLength(graph,path)>MIN_OUTPUT_SEQ)
					{
						FinalPaths_all.put(path,new Pair<Integer>(getSuppCalculation(PathReads.get(path)),0));
						debugMes("the unextended path: "+path+" was added to the final paths, with "+getSuppCalculation(PathReads.get(path)) +" support",15);
					} 
					removePaths.add(path);
				}
			}

			for (List<Integer> path : removePaths)
			{
				debugMes("path "+ path +" wasnt extended and is captured accordingly.",15);
				Paths.get(v).remove(path);
				Extensions.remove(path);
			}
		}
		
		// ***************************
		// end of path constructions
		// ***************************
		

		for (List<Integer> path : Paths.get(T_VERTEX))
		{
			if (getSeqPathLength(graph,path)>MIN_OUTPUT_SEQ)
			{
				FinalPaths_all.put(path,new Pair<Integer>(getSuppCalculation(PathReads.get(path)),0));
				if (path.get(0).intValue() == ROOT.getID())
					debugMes("the finished path: "+ path+" was added to the final paths, with "+getSuppCalculation(PathReads.get(path))+" support",15);
				else
					debugMes("the finished (from middle unextended) path: "+ path+" was added to the final paths, with "+getSuppCalculation(PathReads.get(path)) +" support",15);
			}
		}


		if (FinalPaths_all.size() > 1)
			FinalPaths_all = remove_identical_subseqs(FinalPaths_all, graph, PathReads);
		
	

		return new Pair<HashMap<List<Integer>, Pair<Integer>>>(FinalPaths_diff,FinalPaths_all);

	}



	private static boolean parents_all_visited(SeqVertex v,
			HashMap<Integer, Boolean> node_ID_visited,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		
		for (SeqVertex pred : graph.getPredecessors(v)) {
			if (pred.getID() > 0 && ! node_ID_visited.containsKey(pred.getID())) {
				return(false);
			}
		}
		
		return(true);
		
	}


	private static List<List<Integer>> remove_lesser_supported_paths_EM(
			List<List<Integer>> all_paths,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> pathReads, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			PathExpressionComparator pc, 
			HashMap<List<Integer>, Integer> separate_gene_ids) {
		
		
		
		// determine max expr per gene
		HashMap<Integer,Float> gene_to_max_expr = new HashMap<Integer,Float>();
		
		for (List<Integer> path : all_paths) {
			
			if (separate_gene_ids.containsKey(path)) {

				Integer gene_id = separate_gene_ids.get(path);
				
				float expr = pc.get_expr(path);
				if (gene_to_max_expr.containsKey(gene_id)) {
					if (gene_to_max_expr.get(gene_id) < expr) {
						gene_to_max_expr.put(gene_id, expr);
					}
				}
				else {
					gene_to_max_expr.put(gene_id, expr);
				}
			}
		}

		
		List<List<Integer>> paths_to_keep = new ArrayList<List<Integer>>();
		
		// retain only those genes that have at least x% expr of the dominant isoform.
		for (List<Integer> path : all_paths) {
			
			boolean keep = true;
		
			if (separate_gene_ids.containsKey(path)) {

				Integer gene_id = separate_gene_ids.get(path);
				
				float expr = pc.get_expr(path);
				float max_iso_expr = gene_to_max_expr.get(gene_id);
				
				float pct_iso_expr = expr / max_iso_expr * 100;
				
				if (pct_iso_expr >= MIN_RELATIVE_ISOFORM_EXPRESSION) {

					// keep it.
					debugMes("Keeping isoform: " + path + " as having " + pct_iso_expr + "% dom. iso expr for gene.", 15);


				}
				else {
					keep = false;
					debugMes("*Excluding isoform: " + path + " as having " + pct_iso_expr + "% dom. iso expr for gene.", 15);

				}
				
			}
			
			if (keep) {
				paths_to_keep.add(path);
			}
			
			
		}
		
		return(paths_to_keep);
			
		
		
	}


	private static List<List<Integer>> remove_lesser_supported_paths(
			List<List<Integer>> paths_ending_at_v,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> pathReads) {
		
		// note, should be already sorted by priority from high-to-low
		
		debugMes("\n## Removing lesser-supported paths that end at V", 10);
		
		HashMap<List<Integer>, List<PairPath>> path_to_compatible_reads = new HashMap<List<Integer>, List<PairPath>>();
		
		// compute compatibility and containments
		for (List<Integer> path : paths_ending_at_v) {
			path_to_compatible_reads.put(path, new ArrayList<PairPath>());
			
			for (PairPath pp : pathReads.get(path).keySet()) {
				if (pp.isCompatible(path)) {
					path_to_compatible_reads.get(path).add(pp);
				}
			}
		}
		
		
		
		List<List<Integer>> best_paths = new ArrayList<List<Integer>>();
		
		// examine them hierarchically and see if lower supported paths continue to add any unique read content
		HashSet<PairPath> all_PairPaths = new HashSet<PairPath>();
		for (List<Integer> path : paths_ending_at_v) {
			List<PairPath> compat_reads = path_to_compatible_reads.get(path);
			
			int count_unique = 0;
			for (PairPath pp : compat_reads) {
				if (! all_PairPaths.contains(pp)) {
					count_unique++;
					all_PairPaths.add(pp);
				}
			}
			
			debugMes("Unique contribution of pairpath: " + count_unique + " of total: " + compat_reads.size() + " from path: " + path, 10);
			if (count_unique > 0)
				best_paths.add(path);
			else
				debugMes("\tdiscarding path due to lack of unique read (pairpath) content: " + path, 10);
			
		}
		
		
		
		return(best_paths);
	}


	private static HashMap<List<Integer>, Pair<Integer>> verifyTripletSupportAcrossPaths(
			HashMap<List<Integer>, Pair<Integer>> finalPaths_all,
			HashMap<Integer, List<List<Integer>>> tripletMapper,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		
		
		HashMap<List<Integer>, Pair<Integer>> triplet_reinforced_paths = new HashMap<List<Integer>, Pair<Integer>>();
		
		
		for (List<Integer> path : finalPaths_all.keySet()) {


			if (path.size() < 3) {
				continue;
			}

			debugMes("Verifying triplets for path: " + path, 10);

			// iterate through triplets
			
			ArrayList<Pair<Integer>> triplet_reinforced_regions = new ArrayList<Pair<Integer>>();
			
			int begin = 0;
			
			for (int i = 1; i < path.size()-1; i++) {
				Integer central_id = path.get(i);
				Integer left_id = path.get(i-1);
				Integer right_id = path.get(i+1);

				List<Integer> adjacency_path = new ArrayList<Integer>();
				adjacency_path.add(left_id);
				adjacency_path.add(central_id);
				adjacency_path.add(right_id);


				if (tripletMapper.containsKey(central_id)) {
					debugMes("triplet adjacency_path of node: " + central_id + " => " + adjacency_path + "OK", 10);
				}
				else {
					debugMes("triplet adjacency_path of node: " + central_id + " => " + adjacency_path + "*** MISSING ***", 10);
					triplet_reinforced_regions.add(new Pair<Integer>(begin, i));
					begin = i;
				}
			}
			triplet_reinforced_regions.add(new Pair<Integer>(begin, path.size()-1));

			for (Pair<Integer> subpath_range : triplet_reinforced_regions) {
				
				
				
				Integer start_node_index = subpath_range.getFirst();
				Integer stop_node_index = subpath_range.getSecond();
				
				debugMes("Processing Triplet-Reinforced region: " + path.subList(start_node_index, stop_node_index + 1), 10);
				
				// see if the start node looks like a hub
				Integer start_node_id = path.get(start_node_index);
				SeqVertex start_node = getSeqVertex(graph, path.get(start_node_index));
				if (start_node_id >= 0 && graph.getSuccessorCount(start_node) > 1) {
					start_node_index++;
				}
				
				Integer stop_node_id = path.get(stop_node_index);
				SeqVertex stop_node = getSeqVertex(graph, stop_node_id);
				if (stop_node_id >= 0  && graph.getPredecessorCount(stop_node) > 1) {
					stop_node_index--;
				}
				
				if (start_node_index <= stop_node_index) {
					
					List<Integer> refined_triplet_path = path.subList(start_node_index, stop_node_index+1);
					debugMes("Refined triplet-reinforced path= " + refined_triplet_path, 10);
					debugMes("Start node: " + refined_triplet_path.get(0) + " has successor count: " + graph.getSuccessorCount(getSeqVertex(graph, refined_triplet_path.get(0))), 10);
					debugMes("End node: " + refined_triplet_path.get(refined_triplet_path.size()-1) + " has predecessor count: " + graph.getPredecessorCount(getSeqVertex(graph, refined_triplet_path.get(refined_triplet_path.size()-1))), 10);
					
					refined_triplet_path = ensure_path_has_sinks(refined_triplet_path);
					
					triplet_reinforced_paths.put(refined_triplet_path, new Pair<Integer>(1,1));
					
				}
			}
			
		}
			
		return(triplet_reinforced_paths);
		
	}


	/**
	 * given these paths, and reads, re-calc the FPKM of each path
	 * @param FinalPaths
	 * @param PathReads
	 */
	private static void illustrateFinalPaths(
			HashMap<List<Integer>, Pair<Integer>> FinalPaths,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {


		for (List<Integer> path : FinalPaths.keySet())
		{


			debugMes("\nPATH: " + path, 15);

			Integer supp = 0;
			Integer totalCounts = 0;

			HashMap<PairPath, Integer> containedReads = PathReads.get(path);
			
			String ascii_illustration = getPathMappingAsciiIllustration(path, containedReads);
			debugMes("\nPath Illustration:\n\n" + ascii_illustration + "\n", 5);

			//TODO: enable printing at lower verbose level, but note that crazy long paths cause serious performance problems for generating these illustrations... some fine-tuning definitely required there.



		}

	

	}


	/**
	 * Go over all final paths, and combine those that are too similar.
	 * @param graph
	 * @param FinalPaths
	 * @param PathReads
	 * @param topOrderInts 
	 */
	
	/*
	
	private static void combineSimilarPaths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<List<Integer>, Pair<Integer>> FinalPaths,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {


		List<List<Integer>> removeSimilarPaths = new ArrayList<List<Integer>>();

		debugMes("\n\n===========\nmethod: combineSimilarPaths()", 15);

		Iterator<List<Integer>> i1,i2;
		String path1S="", path2S="";
		int index1;
		int pathcount1 = 0;
		int pathcount2 = 0;

		for (i1=FinalPaths.keySet().iterator() ; i1.hasNext() ; )
		{
			List<Integer> path1 = i1.next();
			path1S = getPathSeq(graph, path1);

			pathcount1++;

			boolean gotToi1 = false;
			for (i2=FinalPaths.keySet().iterator() ; i2.hasNext() ; )
			{

				pathcount2++;

				debugMes("\n************\nComparing paths: (" + pathcount1 + "," + pathcount2 + ")", 15);

				List<Integer> path2 = i2.next();

				while (!gotToi1 && i2.hasNext())
				{
					if (path2.equals(path1))
						gotToi1 = true;
					path2 = i2.next();

				}

				if (path2.equals(path1))
					break;

				index1=path1S.length();

				path2S = getPathSeq(graph, path2);

				boolean noOverlappingVers = true;
				int v1 = -1,v2 = -1, index2=-1;
				for (int j1=path1.size()-1; j1>0 && noOverlappingVers ; j1--)
				{
					v1 = path1.get(j1);
					index1 -= getSeqVertex(graph, v1).getName().length();
					if (v1!=T_VERTEX.getID())
					{
						index2=path2S.length();
						for (int j2=path2.size()-1; j2>0 && noOverlappingVers ; j2--)
						{
							v2 = path2.get(j2);
							index2 -= getSeqVertex(graph, v2).getName().length();

							if (v1==v2)
								//update noOverlappingVers, so we'll get out of the loop:
								noOverlappingVers = false;
						}
					}
				}

				if (!noOverlappingVers) //check only paths that share vertices
				{
					index1 += getSeqVertex(graph, v1).getName().length();
					index2 += getSeqVertex(graph, v2).getName().length();
					debugMes("checking paths: "+ path1+ 
							"(len="+path1S.length()+") and "+path2+"(len="+path2S.length()+")",15);

					if (path1.lastIndexOf(T_VERTEX.getID())==-1)
						index1--;

					if (path2.lastIndexOf(T_VERTEX.getID())==-1)
						index2--;
					if (twoPathsAreTooSimilar(graph, path1, path2))
					{
						debugMes("they are too similar!",15);	
						//remove the shorter path
						removeTheLesserSupportedPath(path1S,path2S,path1,path2,removeSimilarPaths,PathReads);
					}
				}
			}
		}

		for (List<Integer> path2Remove : removeSimilarPaths)
		{
			debugMes("The final path "+path2Remove+" was removed because it was too close to another path",15);
			FinalPaths.remove(path2Remove);
		}
	}

	*/

	/**
	 * check for similar paths that end at V, and start at different nodes
	 * remove the shortest of the two
	 * @param graph
	 * @param v
	 * @param Paths
	 * @param PathReads
	 * @param Extensions
	 * @param topOrderInts 
	 */
	

	
	/*
	
	private static List<List<Integer>> combineSimilarPathsThatEndAtV(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,  SeqVertex v,
			HashMap<SeqVertex, List<List<Integer>>> Paths,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads,
			HashMap<List<Integer>, Boolean> Extensions) {



		// new Throwable().printStackTrace();

		int vertex_id = v.getID();
		int total_num_paths = Paths.get(v).size();

		debugMes("method: combineSimilarPathsThatEndAtV(" + vertex_id + ") with "+total_num_paths+ " paths", 10);
		debugMes("paths are: "+Paths.get(v),17);

		List<List<Integer>> removeSimilarPaths = new ArrayList<List<Integer>>();
		List<Integer> removedPathsIndices = new ArrayList<Integer>();
		String path1S="", path2S="";
		Iterator<List<Integer>> i1, i2;
		int index1, index2, rIndex;

		int pathCount1 = 0;
		int pathCount2 = 0;

		if (total_num_paths<=1)
			return null;


		// all-vs-all comparison among the paths ending at v
		for (i1=Paths.get(v).iterator() ; i1.hasNext() ; )
		{
			List<Integer> path1 = i1.next();
			path1S = getPathSeq(graph, path1);
			index1 = path1S.length()-1;

			pathCount1++;
			pathCount2 = 0;

			if (removedPathsIndices.contains(pathCount1)) {
				continue;
			}

			boolean gotToi1 = false;
			for (i2=Paths.get(v).iterator() ; i2.hasNext() ; )
			{
				List<Integer> path2 = i2.next();

				pathCount2++;

				debugMes("\r*V[" + vertex_id + "] Comparing " + total_num_paths + " paths, pairs:(" + pathCount1 + "," + pathCount2 + ")   ", 15);


				while (!gotToi1 && i2.hasNext())
				{
					if (path2.equals(path1))
						gotToi1 = true;
					path2 = i2.next();
					pathCount2++;
				}

				if (path2.equals(path1))
					break;

				// one of these paths were removed already
				if (removedPathsIndices.contains(pathCount2)) {
					continue;
				}	

				

				path2S = getPathSeq(graph, path2);
				index2 = path2S.length()-1;
				debugMes("checking for similarity the two paths: "+path1+ 
						"(len="+path1S.length()+");"+path2+"(len="+path2S.length()+")",15);

				if (twoPathsAreTooSimilar(graph, path1, path2))
				{
					debugMes("they are too similar!",15);	
					//remove the shorter path
					rIndex = removeTheLesserSupportedPath(path1S,path2S,path1,path2,removeSimilarPaths,PathReads);
					if (rIndex == 1)// the first path was removed
						removedPathsIndices.add(pathCount1);
					else
						removedPathsIndices.add(pathCount2);
				}
			}
		}

		for (List<Integer> path2Remove : removeSimilarPaths)
		{
			debugMes("The path "+path2Remove+" was removed because it was too close to another path",12);

			Paths.get(v).remove(path2Remove);
			Extensions.remove(path2Remove);

		}

		return(removeSimilarPaths);

	}

	
	*/
	

	/** 
	 * compare the sequences of the two paths, and return true if they are more than MIN_PERCENT_IDENTITY_SAME_PATH. 
	 * @param path1s
	 * @param path2s
	 * @param topOrderInts 
	 * @return
	 */
	private static boolean twoPathsAreTooSimilar(
			DirectedSparseGraph<SeqVertex, 
			SimpleEdge> graph,
			List<Integer> path1, 
			List<Integer> path2
			) {

		debugMes("\n\n****\n\nchecking twoPathsAreTooSimilar (" + path1 + "," + path2 + ")\n****\n\n", 15);

		if (! PairPath.haveAnyNodeInCommon(path1, path2)) {
			debugMes("paths: " + path1 + path2 + " have no node in common, cannot be too similar.", 15);
			return(false); // if no node in common, then they shouldn't be too similar.
		}
		
		AlignmentStats numTotalMismatchesAndGaps = getPrevCalcNumMismatches(graph, path1, path2);
		
		
		int shorterLen = Math.min(getSeqPathLength(graph,path1),getSeqPathLength(graph,path2));
		float path_per_id = 100 - (float)numTotalMismatchesAndGaps.mismatches/shorterLen * 100;
		boolean tooSimilar = isThisTooSimilar(numTotalMismatchesAndGaps.mismatches, numTotalMismatchesAndGaps.max_internal_gap_length, path_per_id);

		DecimalFormat df = new DecimalFormat("#.##");

		debugMes("\n\n====\nRunning PATH alignment of : " + path1 + " to " + path2 + " :: numMM:" + numTotalMismatchesAndGaps.mismatches  
				+ ", max_internal_gap: " + numTotalMismatchesAndGaps.max_internal_gap_length
				+ ", path_per_id = " + df.format(path_per_id) + ", tooSimilar: " + tooSimilar, 15);

		debugMes(numTotalMismatchesAndGaps.toString(), 18);

		// compare to doing a full sequence alignment:
		
		if (false) {
			String path1_seq = getPathSeq(graph, path1);
			String path2_seq = getPathSeq(graph, path2);

			Alignment alignment;

			if (SMITH_WATERMAN_ALIGN_FLAG) {

				debugMes("-running Smith-Waterman alignment of path sequences", 15);
				alignment = NWalign.run_SW_alignment("A", path1_seq, "B", path2_seq, 4, -5, 10, 1);
				
			}
			else {
				
				// Needleman Wunsch Global Alignment is default
				debugMes("-running Needleman-Wunsch alignment of path sequences", 15);
				alignment = NWalign.run_NW_alignment("A", path1_seq, "B", path2_seq, 4, -5, 10, 1);   //NW locks up or takes too long with very long sequences (eg. 40kb align to 6kb)


			}
			AlignmentStats a = new AlignmentStats(alignment);
			debugMes("\n\n====\nSEQUENCE_ALIGNMENT_RESULTS:\n" + a.toString(), 15);
			
			debugMes (new jaligner.formats.Pair().format(alignment), 15);
		}
		

		return(tooSimilar);


	}

	/**
	 * for p1, p2 find the latest nodes that they share (v2)
	 * by going backwards on the paths, while using the topological order of the nodes, and advancing while keeping them in order.
	 * @param graph 
	 * @param path1
	 * @param path2
	 * @param topOrderInts 
	 * @return
	 */
	private static Integer findLastSharedNode(DirectedSparseGraph<SeqVertex,SimpleEdge> graph, 
			List<Integer> path1,
			List<Integer> path2) {
		
		path1 = PairPath.trimSinkNodes(path1);
		path2 = PairPath.trimSinkNodes(path2);
	
		if (path1.isEmpty() || path2.isEmpty()) {
			return(-1);
		}
		
		List<SeqVertex> reversePath1 = getReverseSeqVertexPath(graph,path1);
		List<SeqVertex> reversePath2 = getReverseSeqVertexPath(graph,path2);

		Iterator<SeqVertex> p1_iter = reversePath1.iterator();
		Iterator<SeqVertex> p2_iter = reversePath2.iterator();

		SeqVertex p1_v = p1_iter.next();
		SeqVertex p2_v = p2_iter.next();

		SeqVertexFinishTimeComparator finishingTimeComparator = new SeqVertexFinishTimeComparator();


		while (p1_v != p2_v )
		{
			
			if (finishingTimeComparator.compare(p1_v,p2_v)>=0) {
				if (p1_iter.hasNext())
					p1_v = p1_iter.next();
				else 
					break;
			}
			else if (p2_iter.hasNext())
				p2_v = p2_iter.next();
			else
				break;
				
		}
		return (p1_v==p2_v)? p1_v.getID() : -1;
	}



	/**
	 * given the graph and a list of integers, return the reverse list of seqVertices
	 * @param graph
	 * @param path
	 * @return
	 */
	private static List<SeqVertex> getReverseSeqVertexPath(DirectedSparseGraph<SeqVertex,SimpleEdge> graph, List<Integer> path) {
		List<SeqVertex> res = new ArrayList<SeqVertex>();
		for (int i=path.size()-1; i>=0 ; i--){
			res.add(getSeqVertex(graph, path.get(i)));
		}
		return res;
	}






	/** 
	 * given the key of the two paths, return their number of matches.
	 * If this calculation hasn't been done before, calc and save it.
	 * @param key
	 * @return
	 */
	private static AlignmentStats getPrevCalcNumMismatches (			
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			List<Integer> path1, List<Integer> path2) {

		debugMes("getPrevCalcNumMismatches: Path1: " + path1 + "  Path2: " + path2, 15);
		
		
		// Not penalizing end gaps
		boolean is_at_start_of_graph = ( ((! path1.isEmpty()) && path1.get(0) == -1) || ( (! path2.isEmpty()) && path2.get(0) == -1) );
		boolean is_at_end_of_graph =  ( ((! path1.isEmpty()) && path1.get(path1.size()-1) == -2) || ( (! path2.isEmpty()) && (path2.get(path2.size()-1) == -2)) ); 
		
		
		String P1_s = path1+"";
		String P2_s = path2+"";
		Comparator<String> stringComp = String.CASE_INSENSITIVE_ORDER;
		int compRes = stringComp.compare(P1_s, P2_s);
		String key = (compRes>=0)? P1_s+";"+P2_s : P2_s+";"+P1_s;
		
		if (NUM_MISMATCHES_HASH.containsKey(key)) {
			
			AlignmentStats a = NUM_MISMATCHES_HASH.get(key);
			
			// Already computed it, used cached value
			debugMes("key: " + key + ", cached as: " + a.toString(), 15);
			NUM_MISMATCHES_HASH.put(key, a); 
			
			return(a);
		}
		
		
		// both paths are effectively empty
		else if ( (path1.isEmpty() || (path1.size() == 1 && path1.get(0) < 0))
				
				&& 
				
				(path2.isEmpty() || (path2.size() == 1 && path2.get(0) < 0)) 
						) {
			
			
			AlignmentStats a = new AlignmentStats();
		
			NUM_MISMATCHES_HASH.put(key, a); 
			
			return(a);
			
		}
		
		// paths are identical
		else if (path1.equals(path2)) {
			
			AlignmentStats a = new AlignmentStats();
			// perfect matches, no gaps.
			String path1s = getPathSeq(graph, path1);
			a.matches = path1s.length();
			a.alignment_length = path1s.length();
			
			debugMes("paths are equivalent: Path1:" + path1 + ", Path2:" + path2 + " and have alignment stats:" + a.toString(), 15);
			
			NUM_MISMATCHES_HASH.put(key, a);  // cache results
			
			return(a);
			
		}
		
		// empty path 1, but have path 2
		else if (path1.isEmpty() || (path1.size() == 1 && path1.get(0) < 0)) {
			AlignmentStats a = new AlignmentStats();
			
			if (! (is_at_start_of_graph || is_at_end_of_graph)) {
				Integer path2_seq_len = getPathSeq(graph, path2).length();
				a.max_internal_gap_length = path2_seq_len;
				a.gaps = path2_seq_len;
			}
			debugMes("empty path1 vs " + path2 + " = " + a.toString(), 15);
			
			return(a);
			
		}
		
		
		// empty path 2, but have path 1
		else if (path2.isEmpty() || (path2.size() == 1 && path2.get(0) < 0)) {
			AlignmentStats a = new AlignmentStats();
			
			if (! (is_at_start_of_graph || is_at_end_of_graph)) {
				Integer path1_seq_len = getPathSeq(graph, path1).length();
				a.max_internal_gap_length = path1_seq_len;
				a.gaps = path1_seq_len;
			}
			
			debugMes("path1 : " + path1 + " vs empty path2 = " + a.toString(), 15);
			
			return(a);
		}
			
		/*
		else if (path1.get(path1.size() -1) == path2.get(path2.size() -1) ) {
				
			// last elements are the same

			debugMes("paths have same last node: Path1:" + path1 + ", Path2:" + path2, 15);

			AlignmentStats nodeAlignStats = getPrevCalcNumMismatches(graph, 
					path1.subList(path1.size()-1, path1.size()),
					path2.subList(path2.size()-1, path2.size()));

			debugMes("Scores for last node comparison: " + path1 + path2 + nodeAlignStats.toString(), 15);

			// get prefix alignment stats
			List<Integer> subP1_list = path1.subList(0, path1.size()-1);
			List<Integer> subP2_list = path2.subList(0, path2.size()-1);

			AlignmentStats remainingAlignmentStats = getPrevCalcNumMismatches(graph, subP1_list, subP2_list);

			debugMes("prefix alignment stats for: " + subP1_list + subP2_list + remainingAlignmentStats.toString(), 15);

			remainingAlignmentStats  = remainingAlignmentStats.increment_alignment_stats(nodeAlignStats);

			debugMes("summing the alignment scores for : " + path1 + path2 + remainingAlignmentStats.toString(), 15);

			NUM_MISMATCHES_HASH.put(key, remainingAlignmentStats);

			return(remainingAlignmentStats);
		}	
		*/
		
		else {
			
			// No empty path, and alignments haven't been scored/cached yet.

			Integer lastSharedNode = findLastSharedNode(graph,path1,path2);

			if (lastSharedNode!=-1) { // there is a shared node 
				
				
				// run sequence comparison on suffix part, recurse for prefix part if node in common.
				
				debugMes("Paths " + path1 + path2 + " share node " + lastSharedNode, 15);
				
				Integer p1V2index = path1.indexOf(lastSharedNode);
				Integer p2V2index = path2.indexOf(lastSharedNode);

				// get path up to but not including the shared node.
				List<Integer> prefix_path_1 = path1.subList(0, p1V2index);
				List<Integer> prefix_path_2 = path2.subList(0, p2V2index); 

				debugMes("getting prefix alignment for " + prefix_path_1 + prefix_path_2, 15);
				
				// recurse to get the prefix alignment
				AlignmentStats prefixStats = getPrevCalcNumMismatches(graph, prefix_path_1, prefix_path_2);
				debugMes("path prefix alignment stats for: " + prefix_path_1 + " and " + prefix_path_2 + " : " + prefixStats.toString(), 15);

				// get alignment for shared node pair
				List<Integer> shared_node_path1 = path1.subList(p1V2index, p1V2index+1);
				List<Integer> shared_node_path2 = path2.subList(p2V2index, p2V2index+1);
				
				AlignmentStats sharedNodeStats = getPrevCalcNumMismatches(graph, shared_node_path1, shared_node_path2);
				
				// add matches for the current alignment
				List<Integer> suffix_path_1 = new ArrayList<Integer>();
				List<Integer> suffix_path_2 = new ArrayList<Integer>();
				
				if (p1V2index < path1.size()-1) {
					suffix_path_1 = path1.subList(p1V2index+1,path1.size());
				}
				if (p2V2index < path2.size()-1) {
					suffix_path_2 = path2.subList(p2V2index+1,path2.size());
				}
				
				debugMes("getting suffix alignment for: " + suffix_path_1 + suffix_path_2, 15);
				AlignmentStats suffixStats = getPrevCalcNumMismatches(graph, suffix_path_1, suffix_path_2);
				
				debugMes("suffix alignment stats: " + suffixStats, 15);
				
				suffixStats = suffixStats.increment_alignment_stats(sharedNodeStats);
				suffixStats = suffixStats.increment_alignment_stats(prefixStats);
				
				debugMes("combining suffix and prefix alignment stats: " + suffixStats, 15);
				
				NUM_MISMATCHES_HASH.put(key, suffixStats);
				
				return(suffixStats);
				
			}
			else {
				
				// no commonly shared node.
				// do the actual alignment.
				                                        
				debugMes("-no shared node, alignment not cached, computing: " + path1 + " to " + path2, 12);

				String path1s = getPathSeq(graph, path1);
				String path2s = getPathSeq(graph, path2);

				debugMes("-path1s length: " + path1s.length() + ", path2s length: " + path2s.length(), 12);

				//TODO: If one path sequence is a substring of the other, no reason to do an alignment.
				//  this can be known based on the path list comparison, without needing to do a string comparison.

				//align the two seqs

				
				AlignmentStats stats;

				if (   (path1s.length() > MAX_SEQ_LEN_DP_ALIGN && path2s.length() > MAX_SEQ_LEN_DP_ALIGN) 
						||
						path1s.length() > 100000   // problems can arise in the alignment code if either seq is longer
						||
						path2s.length() > 100000
						
						) {
					// heuristic...   if seqs are this long, it's highly unlikely they'd be too similar anyway
					stats = ZipperAlignment.doZipperAlignment("A", path1s, "B", path2s);
				}
				
				// very short sequence in path1 or path2
				else if (  (path1s.length() < 10 && path2s.length() > 20)
							||
							(path1s.length() > 20 && path2s.length() < 10) ) 
				{

					if (is_at_start_of_graph) {
						stats = ZipperAlignment.doZipperAlignmentAnchorRight("A", path1s, "B", path2s);

					}
					else if (is_at_end_of_graph) {
						stats = ZipperAlignment.doZipperAlignmentAnchorLeft("A", path1s, "B", path2s);
					}
					else {
						// let ZipperAlignment try to figure it out
						
						stats = ZipperAlignment.doZipperAlignment("A", path1s, "B", path2s);
					}

				}
				
				
				else if (SMITH_WATERMAN_ALIGN_FLAG) {
					debugMes("-running Smith-Waterman alignment of path sequences", 15);
					Alignment alignment = NWalign.run_SW_alignment("A", path1s, "B", path2s, 4, -5, 10, 1);
					debugMes (new jaligner.formats.Pair().format(alignment), 17);
					stats = new AlignmentStats(alignment);
				}
				else {
					// Needleman Wunsch Global Alignment is default
					debugMes("-running Needleman-Wunsch alignment of path sequences", 15);
					Alignment alignment	= NWalign.run_NW_alignment("A", path1s, "B", path2s, 4, -5, 10, 1);   //NW locks up or takes too long with very long sequences (eg. 40kb align to 6kb)
					debugMes (new jaligner.formats.Pair().format(alignment), 17);
					stats = new AlignmentStats(alignment);
				}


				int alignment_length = stats.alignment_length;
				int matches = stats.matches;
				int mismatches = stats.mismatches;
				int gaps = stats.gaps;

				int right_gap_len = stats.right_gap_length;
				int left_gap_len = stats.left_gap_length;
				int max_internal_gap_length = stats.max_internal_gap_length;

				float percent_A_in_alignment = (float) stats.get_count_of_bases_in_aligned_region("A") / (path1s.length()) * 100;
				float percent_B_in_alignment = (float) stats.get_count_of_bases_in_aligned_region("B") / (path2s.length()) * 100;

				debugMes("Percent A in alignment = " +  stats.get_count_of_bases_in_aligned_region("A") + " / " + path1s.length() + " = " + percent_A_in_alignment + "%",15);
				debugMes("Percent B in alignment = " + stats.get_count_of_bases_in_aligned_region("B") + " / " + path2s.length() + " = " + percent_B_in_alignment + "%",15);

				float max_percent_aligned = Math.max(percent_A_in_alignment, percent_B_in_alignment);


				float percent_identity = (float)matches/(matches+mismatches) * 100;
				float percent_gapped = (float)gaps/alignment_length * 100;

				debugMes("Matches: " + matches + ", Mismatches: " + mismatches + ", gaps: " + gaps + ", align_len: " + alignment_length,15);
				debugMes("percent_identity: " + percent_identity + ", percent_gapped: " + percent_gapped,15);
				debugMes("max_percent_aligned: " + max_percent_aligned,15);
				debugMes("max internal gap length: " + max_internal_gap_length + "\n",15);


				int total_significant_diffs = 0;
				if (is_at_start_of_graph || is_at_end_of_graph) {
					total_significant_diffs = mismatches + gaps; 
					debugMes("(start of graph) Total number of significant alignment diffs = (mismatches: " + mismatches 
							+ " + internal_gaps: " + gaps
							+ " + right_gap_len: "+ right_gap_len
							+ "  = " + total_significant_diffs, 15); 

					// the max internal gap length value based ignores the left gap length
					if (is_at_start_of_graph) {
						stats.left_gap_length = 0; 
						if (! is_at_end_of_graph) {
							// deal with right-gap in alignment stats
							stats.max_internal_gap_length = Math.max(stats.max_internal_gap_length, stats.right_gap_length);
							total_significant_diffs += stats.right_gap_length;
							stats.gaps += stats.right_gap_length;
						}
					}
					if (is_at_end_of_graph) {
						stats.right_gap_length = 0;
						if (! is_at_start_of_graph) {
							// deal with left-gap in alignment stats
							stats.max_internal_gap_length = Math.max(stats.max_internal_gap_length, stats.left_gap_length);
							total_significant_diffs += stats.left_gap_length;
							stats.gaps += stats.left_gap_length;
						}
					}


				}
				else {
					total_significant_diffs = mismatches + gaps + left_gap_len + right_gap_len; // all gaps count TODO: ignore right gap length if at end of graph
					debugMes("(internal of graph) Total number of significant alignment diffs = (mismatches: " + mismatches 
							+ " + internal_gaps: " + gaps
							+ " + left_gap_len: " + left_gap_len 
							+ " + right_gap_len: "+ right_gap_len
							+ "  = " + total_significant_diffs, 15); 


					// adjust max internal gap length value based on left or right gap lengths, since this is an internal node
					stats.max_internal_gap_length = Math.max(stats.max_internal_gap_length, stats.left_gap_length);
					stats.max_internal_gap_length = Math.max(stats.max_internal_gap_length, stats.right_gap_length);

				}

				stats.total_not_matched = total_significant_diffs; // update based on above.

				debugMes("AlignmentStats: " + stats.toString(), 15);
				
				NUM_MISMATCHES_HASH.put(key, stats); 

				return NUM_MISMATCHES_HASH.get(key);
			}
		}

	}

	/**
	 * given all the params, decide if the two seqs are too similar
	 * FIXME - find a better criteria.
	 * @param numMM - number of mismatches
	 * @param longestMMstretch
	 * @param shortestLen
	 * @return
	 */




	private static boolean isThisTooSimilar(int numMM, int max_internal_gap_length, float percent_identity) { // number of differences, not just mismatches, includes gaps


		DecimalFormat df = new DecimalFormat("#.##");



		boolean too_similar = ( max_internal_gap_length <= MAX_INTERNAL_GAP_SAME_PATH
				&&
				(	numMM <= MAX_DIFFS_SAME_PATH || percent_identity >= MIN_PERCENT_IDENTITY_SAME_PATH));


		debugMes("the two paths have these stats: numMM="+numMM
				+ ", max_internal_gap_length=" + max_internal_gap_length
				+  ", identity="+df.format(percent_identity)+"%" 
				+ ", tooSimilar: " + too_similar,15);

		return (too_similar);  // same as saying they are too similar... I just process the logic better in the terms of them not being too different.

	}


	/**
	 * given two paths (and their seqs) remove the shorter path, and add its reads to the other one.
	 * if the are equal in length, remove the lighter one.
	 * @param path1S
	 * @param path2S
	 * @param path1
	 * @param path2
	 * @param removeSimilarPaths
	 * @param PathReads
	 */
	private static int removeTheLesserSupportedPath(String path1S, String path2S,
			List<Integer> path1, List<Integer> path2, List<List<Integer>> removeSimilarPaths,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {
		List<Integer> path2remove,path2keep;

		int sum1=0,sum2=0;
		if (PathReads.get(path1)!=null)
			for (Integer s : PathReads.get(path1).values())
				sum1+=s;
		if (PathReads.get(path2)!=null)
			for (Integer s : PathReads.get(path2).values())
				sum2+=s;

		debugMes("Scoring paths based on read support.  Path: " + path1 + " has " + sum1 + " read support, and " + path2 + " has " + sum2 + " read support.", 15);
		
		
		if (sum1<sum2)
		{
			path2remove = path1;
			path2keep = path2;
		} 
		else if (sum1 > sum2)
		{
			path2remove = path2;
			path2keep = path1;
		}
		else {
			// same read support, so choose the longer sequence over the shorter one.

			if (path1S.length() >= path2S.length())
			{
				path2remove = path2;
				path2keep = path1;
			}
			else
			{
				path2remove = path1;
				path2keep = path2;
			}


		}

		debugMes("removing path "+path2remove+" and keeping path "+path2keep,15);

		
		
		if (!removeSimilarPaths.contains(path2remove))
			removeSimilarPaths.add(path2remove);
		
		if (PathReads.get(path2remove)!=null)
		{
			if (PathReads.get(path2keep)==null)
				PathReads.put(path2keep, new HashMap<PairPath,Integer>());

			
			// no longer assuming ownership of the other's reads, as this causes problems!
			//PathReads.get(path2keep).putAll(PathReads.get(path2remove));
			
			
			PathReads.remove(path2remove);
		}	

		return (path2remove==path1)? 1:2;
	}




	/** 
	 * Given this path, ask whether it has enough support, either by last triplet, or by length
	 * @param readsOfPathUntilV - reads of this path, so far
	 * @param path - the path so far
	 * @param u - the extension to the path
	 * @param graph
	 * @param dijkstraDisWoVer
	 * @return
	 */
	private static boolean pathHasEnoughReadSupport(
			HashMap<PairPath, Integer> readsOfPathUntilV, 
			List<Integer> path,
			SeqVertex u, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer) {


		List<Integer> pathWU = new ArrayList<Integer>(path);
		pathWU.add(u.getID());

		List<Integer> subPath = new ArrayList<Integer>();
		subPath.add(0, u.getID());

		SeqVertex v = getSeqVertex(graph,path.get(path.size()-1));



		if (LENIENT_PATH_CHECKING) {

			// nodes u and v exist within a read pair path, and the read(s) are compatible with the this tentative path.
			return(pathHasTerminalCompatibleReadSupport(path, v, u, graph, readsOfPathUntilV, dijkstraDisWoVer));

		}

		else if (USE_TRIPLETS)  // never do it this way, option turned off permanently but retained for legacy sake.
		{

			subPath.add(0, v.getID());
			if (path.size()>1)
				subPath.add(0,path.get(path.size()-2));



			return (subPathHasEnoughReadSupport(pathWU, readsOfPathUntilV, subPath, graph, dijkstraDisWoVer));
		}

		else{
			// default method

			int lookBack = PATH_REINFORCEMENT_DISTANCE; 
			int lenSoFar = u.getNameKmerAdj().length();
			for (int j = path.size()-1 ; j>=0 && lenSoFar < lookBack; j--){
				SeqVertex vLast = getSeqVertex(graph, path.get(j));
				subPath.add(0, vLast.getID());
				lenSoFar += vLast.getNameKmerAdj().length();
			}

			return (subPathHasEnoughReadSupport(pathWU, readsOfPathUntilV, subPath, graph, dijkstraDisWoVer));
		}
	}



	/**
	 * Check that the given sub-path has N supporting reads or more.
	 * A supporting read is a read that enforces this triplet
	 * @param readsOfPathUntilV
	 * @param subPath
	 * @param graph
	 * @param dijkstraDisWoVer
	 * @return
	 */
	private static boolean subPathHasEnoughReadSupport(
			List<Integer> fullPathWU,
			HashMap<PairPath, Integer> readsOfPathUntilV,
			List<Integer> subPath,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer) {

	
		
		int local_debug_level = 17;
		
		debugMes("-checking if subPath has enough read support. Exploring sub path: " + subPath, local_debug_level);

		// note, subpath contains U as final node.

		
		
		Integer last_subpath_id = subPath.get(subPath.size()-1);
		Integer first_subpath_id = subPath.get(0);
		if (first_subpath_id < 0)
			first_subpath_id = subPath.get(1); // dont want a sink node

		int numberReadsSupporting = 0;
		for (PairPath pPath : readsOfPathUntilV.keySet())
		{

			debugMes("-readsOfPathUntilV: " + pPath, local_debug_level);
			
			boolean thisReadOK = true;

			if (ORIGINAL_PATH_EXTENSIONS) {
				

				// Examining within the context of the entire graph

				for (Integer vTempID : subPath) {
					if (thisReadOK)
						thisReadOK = thisReadOK && 	
						readEnforcesVertex(graph, dijkstraDisWoVer, pPath, getSeqVertex(graph, vTempID));
				}
			}
			else {
				
				
				
				if (! (
						pPath.containsID(last_subpath_id) 
						//&& (first_subpath_id == -1 || pPath.containsID(first_subpath_id) ))  // triplet came in with a sink node attached, not a true triplet //FIXME: should check this earlier than this routine.
						
						&& pPath.isCompatible(fullPathWU)
						
						&& pPath.node_is_contained_or_possibly_in_gap(first_subpath_id, graph, dijkstraDisWoVer)
					  )
					)
				{
					// require pPath contains the first and last ID of the subpath
					thisReadOK = false;
				}
				
				debugMes("\t-checking if pp: " + pPath + " supports extension of " + fullPathWU + " => " + thisReadOK, local_debug_level);
				
				/*
				// COMPATIBLE_PATH_EXTENSIONS MODE, NOW THE DEFAULT
				boolean subPathContained = pPath.containsSubPath(subPath);
				boolean pathWUcompatible = pPath.isCompatible(fullPathWU);

				debugMes("CPATEXT: subPath: " + subPath + " contained by read: " + pPath.get_paths() + " : " + subPathContained, local_debug_level);
				debugMes("CPATEXT: pathWU: " + fullPathWU + " compatible with read: " + pPath.get_paths() + " : " + pathWUcompatible, local_debug_level);

				thisReadOK = (subPathContained && pathWUcompatible);
				 */	
			}
		

			debugMes("examining subPath: " + subPath + " for reinforcement by read: " + pPath.get_paths() + " :" + thisReadOK, local_debug_level);

			if (thisReadOK)
			{
				numberReadsSupporting+=readsOfPathUntilV.get(pPath);
				debugMes("the read "+pPath+"("+readsOfPathUntilV.get(pPath)+") enforces the sub-path ("+subPath+")",local_debug_level);
				if (numberReadsSupporting >= MIN_READ_SUPPORT_THR) {
					break; // no reason to spend time searching more.
				}
			} else
				debugMes("the read "+pPath+"("+readsOfPathUntilV.get(pPath)+") does not enforce the sub-path ("+subPath+")",local_debug_level);

		}

		debugMes("-found: " + numberReadsSupporting + " reads supporting subpath.", local_debug_level);

		boolean res = (numberReadsSupporting>=MIN_READ_SUPPORT_THR);
		if (res)
			debugMes("the sub-path ("+subPath+") has PASSED", local_debug_level);
		else
			debugMes("the sub-path ("+subPath+") has NOT PASSED",local_debug_level);

		return res;	
	}






	/**
	 * Check whether there are at least N reads enforcing 
	 * @param graph
	 * @param dijkstraDis
	 * @param pPath
	 * @param v
	 * @return
	 */
	private static boolean readEnforcesVertex(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer,
			PairPath pPath, SeqVertex v) {

		// read contains vertex v
		if (v==null || pPath.containsID(v.getID()) || v.equals(ROOT) || v.equals(T_VERTEX))
			return true;

		// is path from the root (sink) to first vertex of the read disabled by removing vertex v of subpath?
		SeqVertex firstV = getSeqVertex(graph, pPath.getFirstID());
		if (dijkstraDisWoVer.getDistanceWoVer(ROOT, firstV,v)==null)
			return true;
		// is path from last vertex of the read to the terminal node (sink) disabled by removing vertex v of subpath?
		SeqVertex lastV = getSeqVertex(graph, pPath.getLastID());
		if (dijkstraDisWoVer.getDistanceWoVer(lastV, T_VERTEX,v)==null)
			return true;

		if (pPath.hasSecondPath())
		{

			// is path from beginning to end of read's pair-path disrupted by removing vertex v of subpath?

			//last of first path
			lastV = getSeqVertex(graph, pPath.getLastID_path1());
			//first of second path
			firstV = getSeqVertex(graph, pPath.getFirstID_path2());
			if (dijkstraDisWoVer.getDistanceWoVer(lastV,firstV,v)==null)
				return true;
		}
		return false;
	}


	private static boolean pathHasTerminalCompatibleReadSupport(
			List<Integer> path,
			SeqVertex v, SeqVertex u, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<PairPath, Integer> readsOfPathUntilV,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer) {


		List<Integer> tentativePath = new Vector<Integer>(path);
		tentativePath.add(u.getID());

		Integer v_id = v.getID();
		Integer u_id = u.getID();

		List<Integer> subPath = new ArrayList<Integer>();
		subPath.add(v_id);
		subPath.add(u_id);

		int num_compatible_paths = 0;

		for (PairPath pPath : readsOfPathUntilV.keySet()) {

			//if (pPath.containsID(v_id) && pPath.containsID(u_id)) {
			if (pPath.containsSubPath(subPath)) {
				debugMes("Checking for compatibility.  Path: " + tentativePath +  " with " + pPath, 18);
				// got both terminal path vertices.  Check for read compatibility.
				if (pPath.isCompatible(path)) {
					debugMes("\tPaths ARE compatible.", 18);
					num_compatible_paths += readsOfPathUntilV.get(pPath);
				}

			}

		}

		debugMes("\t" + num_compatible_paths + " read (pair) paths were found to be compatible.", 18);

		if (num_compatible_paths >= MIN_READ_SUPPORT_THR) { // note, not using this as triplet support here. 
			//TODO: rename triplet support var
			return(true);
		}
		else {
			return(false);
		}

	}



	private static boolean vertexPairHasDiscontinuousPathSupport(
			List<Integer> path,
			SeqVertex v, SeqVertex u, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<PairPath, Integer> readsOfPathUntilV,
			DijkstraDistanceWoVer<SeqVertex, SimpleEdge> dijkstraDisWoVer) {


		debugMes("\n\nCurrent path being checked for LENIENT path extension: " + path, 18);

		debugMes("Performing LENIENT path checking between (v,u):\nv: " + v + "\nu: " + u, 18);

		// look for u-v where u is last node of one pairpath, and v-u is the start of another path
		// find a pairpath that ends in v
		// find another pairpath that starts with v-u

		boolean last_vertex_found_as_v = false;
		boolean first_vertices_found_as_vu = false;

		for (PairPath pPath : readsOfPathUntilV.keySet()) {

			debugMes("\t-pairPath: " + pPath, 18);

			SeqVertex last_vertex = getSeqVertex(graph, pPath.getLastID());
			debugMes("\t-Last vertex: " + last_vertex.getID(), 18);

			if (last_vertex.equals(v)) {
				last_vertex_found_as_v = true;
				debugMes("\t\t-found last vertex as (v)", 18);
			}

			List<Integer> first_path = pPath.getPath1();
			if (first_path.size() > 1) {
				SeqVertex first_vertex = getSeqVertex(graph, first_path.get(0));
				SeqVertex second_vertex = getSeqVertex(graph, first_path.get(1));

				debugMes("\t-First,Second: " + first_vertex.getID() + "," + second_vertex.getID(), 18);

				if (first_vertex.equals(v) && second_vertex.equals(u)) {
					first_vertices_found_as_vu = true;
					debugMes("\t\t-found first vertices as (vu)", 18);
				}

			}


			if (first_vertices_found_as_vu && last_vertex_found_as_v) {
				debugMes("\t* FOUND LENIENT EXTENSION", 18);
				return(true);
			}
		}

		debugMes("\t* no LENIENT extension possible", 18);
		return(false);  // no evidence for discontinous support.

	}


	/**
	 * Check whether the pairPath is consistent with the node i
	 * @param pPath
	 * @param i
	 * @param graph
	 * @param dijkstraDis
	 * @return
	 */
	private static boolean readIsConsistentWithNode(PairPath pPath, Integer i,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis) {

		//		if (isReadCircular(graph, pPath))
		//			return false;

		if (pPath.containsID(i) || i<0)
			return true;
		
		SeqVertex vI = getSeqVertex(graph, i);
		SeqVertex firstV = getSeqVertex(graph, pPath.getFirstID());

		// i --> firstV
		if (SeqVertex.isAncestral(vI, firstV, dijkstraDis)>0)
			return true;

		SeqVertex lastV = getSeqVertex(graph, pPath.getLastID());

		// lastV --> i
		if (SeqVertex.isAncestral(lastV,vI,dijkstraDis)>0)
			return true;

		if (pPath.hasSecondPath())
		{
			
			// see if node could be internal to the pair path
			
			//last of first path
			lastV = getSeqVertex(graph, pPath.getLastID_path1());
			//first of second path
			firstV = getSeqVertex(graph, pPath.getFirstID_path2());

			// lastV --> i --> firstV
			if (SeqVertex.isAncestral(lastV,vI,dijkstraDis)>0 && SeqVertex.isAncestral(vI, firstV, dijkstraDis)>0)
				return true;
		}
		
		// not compatible if got here.
		return false;
	}

	
	
	

	/**
	 * given the graph and a list of nodes, calc the length of the seq of this path
	 * @param graph
	 * @param path
	 * @return
	 */
	private static int getSeqPathLength(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, List<Integer> path) {
		
		/*
		int len = 0;
		for (Integer vid : path)
			if (vid>=0)
				len +=getSeqVertex(graph, vid).getNameKmerAdj().length();	
		return len;
		*/
		
		String pathSeq = getPathSeq(graph, path);
		return(pathSeq.length());
		
	}

	/**
	 * return the number of paths
	 * @param paths
	 * @return number of paths
	 */
	private static int getPathsSize(
			HashMap<SeqVertex, List<List<Integer>>> paths) {
		int res = 0;
		for (SeqVertex key : paths.keySet())
		{
			res+=paths.get(key).size();
		}
		return res;
	}


	/**
	 * returns true iff these two nucleotides are equal
	 * @param n1
	 * @param n2
	 * @return
	 */
	private static boolean areTwoNucleotidesEqual(String n1, String n2)
	{
		if (n1.equals(n2))
			return true;

		if (USE_DEGENERATE_CODE && 
				((DEGENERATE_CODE_REV.containsKey(n1) && DEGENERATE_CODE_REV.get(n1).contains(n2)) ||
						(DEGENERATE_CODE_REV.containsKey(n2) && DEGENERATE_CODE_REV.get(n2).contains(n1))))
			return true;

		return false;
	}

	/**
	 * return the degenerate code representation of the given key
	 * @param key
	 * @return
	 * @throws Exception
	 */
	private static String getDegenerateRepresentation(String key) throws Exception {
		if (DEGENERATE_CODE.containsKey(key))
			return DEGENERATE_CODE.get(key);
		else
			throw new Exception("the letters "+key+" do not have a degenerate representation\n");
	}



	/**
	 * sum the counts of all the reads in this hash
	 * @param readHash
	 * @return
	 */
	private static Integer getSuppCalculation(HashMap<PairPath, Integer> readHash) {
		Integer res = 0;
		for (PairPath key : readHash.keySet())
			res = res + readHash.get(key);
		return res;
	}



	/**
	 * Given the new path (with u), and the set of reads that supported the path until v
	 * update the set of reads that support the new path
	 * @param PathReads 
	 * @param pathContainedReads 
	 * @param pathWu
	 * @param readsOfPathUntilV
	 * @param i 
	 * @param dijkstraDis 
	 * @param graph 
	 */
	private static void updateReadsOfPath(HashMap<List<Integer>,HashMap<PairPath,Integer>> PathReads, 
			HashMap<List<Integer>, HashSet<PairPath>> PathContainedReads, List<Integer> pathWu,
			HashMap<PairPath, Integer> readsOfPathUntilV, Integer i, DirectedSparseGraph<SeqVertex, SimpleEdge> graph, DijkstraDistance<SeqVertex, SimpleEdge> dijkstraDis) {

		debugMes("updateReadsOfPath: " + pathWu, 17);
		
		List<Integer> pathMinusU = new ArrayList<Integer>(pathWu);
		pathMinusU.remove(pathMinusU.size()-1); // remove the U o fthe pathWu
		
		// init read path list as needed.
		if (!PathReads.containsKey(pathWu))
			PathReads.put(pathWu, new HashMap<PairPath,Integer>());

		if (! PathContainedReads.containsKey(pathWu))  {
			PathContainedReads.put(pathWu, new HashSet<PairPath>());
		}
	   
		
		int count_total = 0;
		int count_contained_propagated = 0;

		for (PairPath pPath : readsOfPathUntilV.keySet())
		{
			count_total++;
			
			if (!PathReads.get(pathWu).containsKey(pPath)) {  // only if this read doesn't exist in the PathReads for this pathWu
				// if this read is consistent with pathWu, then add it
				//if (readIsConsistentWithNode(pPath,i,graph,dijkstraDis))

				if (PathContainedReads.get(pathMinusU).contains(pPath)) {
					// then PathWU must contain it as well.
					PathContainedReads.get(pathWu).add(pPath);
					PathReads.get(pathWu).put(pPath,readsOfPathUntilV.get(pPath));
					//debugMes("path is contained by pathMinusU: " + pPath, 10);
					
					count_contained_propagated++;
				}

				else {
					if (pPath.isCompatible(pathWu)) { // semi-expensive operation
						if (pPath.isCompatibleAndContainedBySinglePath(pathWu)) {
							PathContainedReads.get(pathWu).add(pPath);
						}
						debugMes("read "+pPath+" is consistent with "+i, 17);
						PathReads.get(pathWu).put(pPath,readsOfPathUntilV.get(pPath));
					}else{
						debugMes("read "+pPath+" is not consistent with "+i, 17);
					}
				}
			}
		}

		float pct_contained_propagated = (float) count_contained_propagated/count_total*100;
		
		debugMes("pct_contained_propagated: " + pct_contained_propagated + "%", 10);
		
	}

	//	/**
	//	 * return true iff this read is circular. 
	//	 * A read is considered circular if its gap includes a circle 
	//	 * (the vertex at the end of path1 is inside a circle or the first vertex of path 2 is inside a circle).
	//	 * @param graph
	//	 * @param readPath
	//	 * @return
	//	 */
	//	private static boolean isReadCircular(DirectedSparseGraph<SeqVertex, SimpleEdge> graph, PairPath readPath)
	//	{
	//		if (!readPath.hasSecondPath())
	//			return false;
	//		
	//		if (getSeqVertex(graph, readPath.getLastID()).isInCircle() ||        //lastID of first path is circular 
	//				getSeqVertex(graph, readPath.getFirstID_path2()).isInCircle()) // firstID of second path is circular
	//				{
	//					debugMes("the read "+readPath+" is circular",10);
	//					return true;
	//				} else
	//					return false;
	//	}


	/**
	 * Return the reads, hashed by their starting vertex
	 * @param graph
	 * @param filename
	 * @param originalVerIDsMapping
	 * @param rootIDs
	 * @param originalGraphKmerToNodeID 
	 * @return
	 * @throws IOException
	 */
	private static HashMap<String, List<Read>> getReadStarts(DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			String filename,
			HashMap<Integer, LocInGraph> originalVerIDsMapping, 
			Vector<Integer> rootIDs, 
			HashMap<String, Integer> originalGraphKmerToNodeID) 
					throws IOException {
		
		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 

		HashMap<String, List<Read>> readNameHash = new HashMap<String, List<Read>>();
		String l = fileB.readLine(); // read header of component
		int numReadsNotMapped = 0;
		int numReadsMapped = 0;
		int line_counter = 0;

		while (fileB.ready())

		{
			l = fileB.readLine();
			if (l.isEmpty())
				continue;

			line_counter++;
			if (line_counter % 1000 == 0 && BFLY_GLOBALS.VERBOSE_LEVEL >= 10)
				System.err.print("\rmapped read [" + line_counter + "]");	
			else if (BFLY_GLOBALS.VERBOSE_LEVEL >= 11) {
				System.err.print("\rmapped read [" + line_counter + "]");
			}

			//			Component 0
			//			>@42MRYAAXX100104:7:100:1000:103#0      11      101393  36      101418          GAAAGACTGTCACCCTTGAGGTGGAGTCCTCTGACACTATTGACAATGTCAAGAGCAAAATCCAAGACAAGGAAGG
			debugMes("Read: " + l, 20);
			String[] fields = l.split("\t");
			fields[0] = fields[0].replaceFirst(">", "");

			List<Integer> pathIDS = null;
			Read r = new Read();
			pathIDS = readAndMapSingleRead(fields,originalVerIDsMapping,graph,r,false,originalGraphKmerToNodeID);

			//debugMes("Threaded Read As: " + r.getName() + " : " + pathIDS, 19);
			
			if (pathIDS==null || (pathIDS!=null && pathIDS.isEmpty()))
			{
				numReadsNotMapped++;
				debugMes("Read could not be threaded: " + r.getName(), 12);
			}else
			{

				//add to readNameHash
				if (!readNameHash.containsKey(r.getName()))
					readNameHash.put(r.getName(), new ArrayList<Read>());

				readNameHash.get(r.getName()).add(r);
				numReadsMapped++;

				//System.err.println(r.getName());
				debugMes("Threaded Read as: " + r.getName() + " : " + pathIDS, 17);
				debugMes("ReadPath@Init: " + r.getName() + " : " + pathIDS, 12);

			}
		}	

		

		//		debugMes("number of reads not found in graph = "+numReadsNotMapped +" of a total of "+(numReadsNotMapped+numReadsMapped),10);
		debugMes("number of reads threaded = "+numReadsMapped 
				+" (from total of "+(numReadsNotMapped+numReadsMapped)+") which came from "
				+ readNameHash.keySet().size() + " pairs",10);

		if (numReadsNotMapped > .5*(numReadsNotMapped+numReadsMapped))
			debugMes("PROBLEM: less than half of the reads were mapped to this graph ("
					+numReadsMapped+"/"+(numReadsNotMapped+numReadsMapped)+")",10);

		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 18) {
			for (String readName : readNameHash.keySet()) {
				String descr = "Read name to pairing info: " + readName + " => "; 
				List<Read> read_list = readNameHash.get(readName);
				for (Read r : read_list) {
					descr += r.getPathIDs();
				}
				debugMes(descr, 15);
			}

		}


		return readNameHash;
	}


	/**
	 * given this read, try and map it to the graph. if rev= true, do it in reverse.
	 * @param fields
	 * @param originalVerIDsMapping
	 * @param graph
	 * @param r 
	 * @param rev
	 * @param originalGraphKmerToNodeID 
	 * @return
	 */
	private static List<Integer> readAndMapSingleRead(String[] fields,
			HashMap<Integer, LocInGraph> originalVerIDsMapping,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, Read r, boolean rev,
			HashMap<String, Integer> originalGraphKmerToNodeID) {

		List<Integer> pathIDS = new ArrayList<Integer>();
		LocInGraph fromV;
		Integer startInRead,endInRead,fromOrigV;

		String name;
		String seq;

		name = fields[0];
		
		if (! TREAT_PAIRS_AS_SINGLE) {
			if (name.endsWith("/1") || name.endsWith("/2") 
					|| name.endsWith("\1") || name.endsWith("\2")
					|| name.endsWith(":1") || name.endsWith(":2")
					)
				name = name.substring(0, name.length()-2);
		}
		
		startInRead = Integer.parseInt(fields[1]);
		//endInRead = Integer.parseInt(fields[3])+KMER_SIZE-1; 
		endInRead = Integer.parseInt(fields[3])+KMER_SIZE; //FIXME: should be as above, but chrysalis appears to be consistently off by one here...
		
		fromOrigV = Integer.parseInt(fields[2]);
		fromV = originalVerIDsMapping.get(fromOrigV);
		seq = fields[6]; //there is an empty field before the seq.

		r.init(name,seq, fromV, startInRead, endInRead, pathIDS);
		

		if (endInRead >= seq.length()) {
			debugMes("read " + name + " has sequence length that is shorter than supposed endInRead marking(" + endInRead + "): " + seq, 0);
			return pathIDS;
		}

		debugMes("Read: " + name + " has start: " + startInRead + ", end: " + endInRead + " and sequence: " + seq, 20);
		
		
		seq = seq.substring(startInRead, endInRead+1);

		debugMes("after extracting substring: " + seq, 20);
		
		
		// in case original node ID was trimmed from graph, try anchoring the sequence from the first
		// recognizable retained node ID
		if (fromV == null) {
			debugMes("Original node ID : " + fromOrigV + " no longer exists ... walking the sequence to try to anchor it to the refined graph:", 20);
			for (int i = 1 + 1; i <= seq.length() - KMER_SIZE; i++) {
				String kmer = seq.substring(i, i+KMER_SIZE);
				if (originalGraphKmerToNodeID.containsKey(kmer)) {
					int ID = originalGraphKmerToNodeID.get(kmer);
					fromV = originalVerIDsMapping.get(ID);
					if (fromV != null) {
						debugMes("Anchored read to graph at position " + (i + startInRead) + " with kmer " + kmer, 20);
						seq = seq.substring(i);
						break;
					}
							
				}
			}
			if (fromV != null) {
				debugMes("recovered mapping of read " + name, 20);
			}
			else {
				debugMes("couldn't recover mapping of read: " + name, 20);
			}
		}
		
		
		
		if (fromV!=null)// && toV!=null)
		{
			Path_n_MM_count best_path = findPathInGraph(graph,seq,fromV,name);
				
			if (best_path != null) {
				
				pathIDS = best_path.path;
				//System.err.println("read path: " + pathIDS);
				
				if (READ_END_PATH_TRIM_LENGTH > 0) {
					// do some read path trimming at beginning and end of path if little support

					pathIDS = best_path.get_trimmed_path(READ_END_PATH_TRIM_LENGTH);
				}
				if (! pathIDS.isEmpty()) {
					r.init(name,seq, fromV, startInRead, endInRead, pathIDS);
					String decorator = (pathIDS.size() < best_path.path.size()) ? " ****** " : "";
					debugMes("Read " + name + " seq " + seq + " threaded as: " + best_path.toString() + ", trimmed to: " + pathIDS + decorator, 17);
				}
				else {
					debugMes("Trimmed path for read: " + name + " threaded as: " + best_path.toString() + " is empty", 19);
					
				}
				
				String pathSeq = getPathSeq(graph, pathIDS);
				
				if (false) { // for debugging 
					if (pathSeq.indexOf(seq) < 0) {
						throw new RuntimeException("Error, read seq: " + seq + "\n does not match threaded seq:\n" + pathSeq);
					}
					else {
						debugMes("Read seq with len=" + seq.length() + " : " + seq 
								+ " matches to path seq with len=" + pathSeq.length() + " : " + pathSeq, 15);
					}
				}

			}

			
		}else
			debugMes("read "+name+" was not mapped to graph. original node doesn't exist anymore ("+fromOrigV+")",19);

		return pathIDS;
	}


	/**
	 * Given the graph, and the read, find the path of the read in the graph
	 * @param graph
	 * @param seq
	 * @param fromV
	 * @return
	 */
	private static Path_n_MM_count findPathInGraph(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, String seq,
			LocInGraph fromV, String readName) {

		
		int local_verbose_level = BFLY_GLOBALS.VERBOSE_LEVEL;
		
		if (readName.startsWith("LR$") && BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
			BFLY_GLOBALS.VERBOSE_LEVEL = 20;
		}
		
		List<Integer> path = new ArrayList<Integer>();

		SeqVertex fromVer = getSeqVertex(graph, fromV.getNodeID());
		
		debugMes("findPathInGraph: V:" + fromV.getNodeID() + ", with seq: " + fromVer.getName(), 20);
		debugMes("trying to start the mapping to node "+fromVer.getID() + " at position: " + fromV.getIndexInNode(), 20);
			
		MAX_MM_ALLOWED_CAP = (int) Math.ceil(seq.length() * MAX_READ_SEQ_DIVERGENCE);
		
		MAX_MM_ALLOWED = MAX_MM_ALLOWED_CAP;
		

		debugMes("\n\nThreading read: " + readName + ", length: " + seq.length()
				+ ", allowing for " + MAX_MM_ALLOWED + " max mismatches.", 17);
		debugMes("Read: " + readName + " sequence is:\n" + seq, 20);
		Integer totalNumMM = 0;
		
		HashMap<String,Path_n_MM_count> best_path_memoization = new HashMap<String,Path_n_MM_count> (); // use DP
		
		Path_n_MM_count best_path_mapping = updatePathRecursively(graph,fromVer.getID(),seq,0, fromV.getIndexInNode(),
												totalNumMM, readName, best_path_memoization);

		
		BFLY_GLOBALS.VERBOSE_LEVEL = local_verbose_level;
		
		if (best_path_mapping != null) {
			
			debugMes("FINAL BEST PATH for " + readName + " is " + best_path_mapping.path + " with total mm: " + best_path_mapping.mismatch_count, 15);
			
			return(best_path_mapping);
		}
		else {
			debugMes("NO_READ_MAPPING_FOUND_FOR: " + readName + "\n\n", 15);
		
			return(null); // no such path found.	
		}
	}
	
	/**
	 * Update the given path recursively
	 * @param path
	 * @param graph
	 * @param fromVers
	 * @param seq
	 * @param locInNode
	 * @param totalNumMM
	 * @param readName
	 * @param best_path_memoization 
	 */
	
	
	private static Path_n_MM_count updatePathRecursively(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Integer fromV_id, String seq, int locInSeq, int locInNode, 
			Integer totalNumMM,String readName, HashMap<String, Path_n_MM_count> best_path_memoization) {

		
		int MIN_SEQ_LENGTH_TEST_DIVERGENCE = 20;
		int MAX_LEFT_END_GAPS = 5;
		
		

		SeqVertex fromV = getSeqVertex(graph, fromV_id);

		Integer numMM = totalNumMM; // init for each node check

		String verSeq = fromV.getName(); //important - full name, not kmer-adjusted name.
		
		
		debugMes("updatePathRecursively(readName=" + readName + 
				", locInSeq: " + locInSeq + " / " + (seq.length() -1) +
				", locInNode: " + locInNode + " / " + (verSeq.length() -1) +
				", totalNumMm: " + totalNumMM, 20);

		
		int startI = locInNode;
		int j=locInSeq, i = startI;

		String read_vertex_start_pos_token = "" + fromV.getID() + "_" + locInNode + "_" + locInSeq;
		
		if (best_path_memoization.containsKey(read_vertex_start_pos_token))  {
			
			Path_n_MM_count best_path = best_path_memoization.get(read_vertex_start_pos_token);
			if (best_path == null) {
				debugMes("MEMOIZATION: indicates this path was a dead end. Not trying again.", 20);
				return(null);
			}
			else {
				debugMes("MEMOIZATION: already stored best path at: " + read_vertex_start_pos_token + " = " + best_path, 20);

				// return a copy, critically important!!! 
				Path_n_MM_count best_path_copy = new Path_n_MM_count(best_path);

				return(best_path_copy);
			}
		}
		debugMes("\ntrying to continue the mapping to node "+ fromV.getShortSeqWID(), 19);

		int length_to_align = Math.min(verSeq.length() - i, seq.length() - j);
		debugMes("-ALIGNING READ SEQ (" + readName + ")\n" + seq.substring(j, j+length_to_align) + " " + j
				   + "\nTo VERTEX (" + fromV.getShortSeqWID() + ") SEQ:\n" + verSeq.substring(i, i+length_to_align) + " " + i, 20);
	
		
		debugMes("Note, rest of read sequence is:\n" + seq.substring(j), 21);

		
		
		
		// zipper align
		boolean failed_alignment = false;
		Integer mm_encountered_here = 0;
		
		for (; i>=0 && i<verSeq.length() && j<seq.length() ; i++,j++)
		{

			String readLetter = ""+seq.charAt(j);
			String verLetter = ""+verSeq.charAt(i);

			String mismatchFlag = (areTwoNucleotidesEqual(readLetter,verLetter)) ? "" : "XXX mismatch XXX";
			debugMes("Comparing read bases: " + i + ":" + readLetter + ", " + j + ":" + verLetter + " " + mismatchFlag, 21);

		
			
			if (!areTwoNucleotidesEqual(readLetter,verLetter)) 
			{
				//we have a mismatch
				numMM++;
				mm_encountered_here++;
				
				if ( (numMM > MAX_MM_ALLOWED)
						||  
						(i >= MIN_SEQ_LENGTH_TEST_DIVERGENCE && (mm_encountered_here/(float)(i)) > MAX_READ_LOCAL_SEQ_DIVERGENCE) 
						)
				{
					failed_alignment = true;
					debugMes("shortcircuiting the zipper test, too many MM or execeeding local seq divergence", 20);
					break; // no point in looking further.
				}

				//TODO: look at mismatch density here as well.

			}
		} // end of mapping read within node

		if (! failed_alignment) {
			
			debugMes("zipper alignment mm: " + mm_encountered_here, 20);
		}
		
		
		
		
		// retain zipper info in case it's better than any DP alignment score
		int zipper_i = i;
		int zipper_j = j;
		int zipper_mm = mm_encountered_here;
		
		
		
		
		// use DP alignment if variation is encountered above. (trying simplest/fastest strategy first)
		
		boolean short_DP_test_passes = true;
		int MIN_LENGTH_TEST_DP = 100;
		
		if (USE_DP_READ_TO_VERTEX_ALIGN && length_to_align > MIN_LENGTH_TEST_DP && mm_encountered_here > 1) {
			
			debugMes("Running short DP test", 20);
			
			j=locInSeq;
			i = startI;
			
			Alignment alignment	= NWalign.run_NW_alignment(		 
					"Vertex", verSeq.substring(i, i+MIN_LENGTH_TEST_DP), 
					"Read", seq.substring(j, j+MIN_LENGTH_TEST_DP), 
					4, -5, 10, 1); 
			
					
			debugMes ("DP test:\n" + new jaligner.formats.Pair().format(alignment), 17);
			AlignmentStats stats = new AlignmentStats(alignment);
			
			mm_encountered_here = stats.mismatches + stats.gaps + stats.left_gap_length;
			
			
			float pct_divergence = mm_encountered_here/(float)(MIN_LENGTH_TEST_DP);
			if ( pct_divergence > MAX_READ_LOCAL_SEQ_DIVERGENCE) {
				debugMes("DP test indicates excessive divergence: " + pct_divergence, 20);
				
				short_DP_test_passes = false;
				
				// leave failed alignment status as is.
			}
			
			// retain earlier zipper stats, regardless of whether or not we go into full DP below.
			i = zipper_i;
			j = zipper_j;
			mm_encountered_here = zipper_mm;
			
		}
		
		int vertex_num_right_end_gaps = 0;
		int read_num_right_end_gaps = 0;
		int max_left_gaps = 0;
		
		
		if (USE_DP_READ_TO_VERTEX_ALIGN && verSeq.length() > 2 && mm_encountered_here > 1 && short_DP_test_passes) {
			
			debugMes("  *Trying again using full DP alignment:", 20);
			
			// reset i and j
			j=locInSeq;
			i = startI;
			
			// try aligning the full vertex sequence w/ extended ref sequence in case it contains small deletions.
			
			
			// Needleman Wunsch Global Alignment is default
			debugMes("-running Needleman-Wunsch alignment of vertex to read", 17);
			/*
			  Alignment alignment	= NWalign.run_NW_alignment(
			 
					"Vertex", verSeq.substring(i, i+length_to_align), 
					"Read", seq.substring(j, j+length_to_align), 
					4, -5, 10, 1);   //NW locks up or takes too long with very long sequences (eg. 40kb align to 6kb)
		    */
			
			int read_length_to_align =  (int) (verSeq.length() * 1.05f);
			if (read_length_to_align + j > seq.length()) {
				read_length_to_align = seq.length() - j;
			}
			
			int bandwidth = (int) (MAX_READ_LOCAL_SEQ_DIVERGENCE * read_length_to_align);
			Alignment alignment	= NWalign.run_NW_banded_alignment(		 
					"Vertex", verSeq.substring(i), 
					"Read", seq.substring(j, j+read_length_to_align), 
					4, -5, 10, 1, bandwidth); 
			
					
			debugMes (new jaligner.formats.Pair().format(alignment), 17);
			AlignmentStats stats = new AlignmentStats(alignment);
			
			mm_encountered_here = stats.mismatches + stats.gaps + stats.left_gap_length + stats.right_gap_length;
			

			
			// check for right end gap in sequence
			String name1 = alignment.getName1();
			char[] vertex_align = alignment.getSequence1();
			char[] read_align = alignment.getSequence2();
			
			if (name1.equals("Read")) {
				char[] swap = vertex_align;
				vertex_align = read_align;
				read_align = swap;
			}
			
			
			
			vertex_num_right_end_gaps = AlignmentStats.get_num_right_end_gaps(vertex_align);
			read_num_right_end_gaps = AlignmentStats.get_num_right_end_gaps(read_align);
			
			max_left_gaps = Math.max(AlignmentStats.get_num_left_end_gaps(vertex_align), 
									 AlignmentStats.get_num_left_end_gaps(read_align));
			
			debugMes("vertex end gaps: " + vertex_num_right_end_gaps, 20);
			debugMes("read end gaps: " + read_num_right_end_gaps, 20);
			
			i = verSeq.length(); // aligning to whole vertex sequence.
			j += read_length_to_align;
			
			if (vertex_num_right_end_gaps > 0) {
				// read extends beyond vertex sequence
				
				
				/*
				 * Vertex          6551 AATTTTCCAGGAAGCTTCAACTGTGTTTGCAAAACTGGATACACAGG---   6597
	                                    |||||||||||||||||||| |||||||||||||||||||||||||.   
	               Read            6430 AATTTTCCAGGAAGCTTCAA-TGTGTTTGCAAAACTGGATACACAGTGTC   6478

				 *  Example of vertex end-gaps.
				 * 
				 */
			
				 j -= vertex_num_right_end_gaps;
			
				mm_encountered_here -= vertex_num_right_end_gaps;
			}
			else if (read_num_right_end_gaps > 0) {
				// vertex extends beyond end of read.
				
				/*
				 * Vertex          2301 TGGGGAAGCAGAACAGTATGTGTGAAGTTTATGTACTGGCACTATAAAAT   2350
                                        |||||||||||||||||||||||||                         
                   Read            2245 TGGGGAAGCAGAACAGTATGTGTGA-------------------------   2269

				 *    Example of read end-gaps
				 */
				
				
				mm_encountered_here -= read_num_right_end_gaps;
			
			}
				
				
			/*
			// if read end gaps, extend to the end of the vertex
			if (read_num_right_end_gaps > 0) {
				
				for (int r = read_num_right_end_gaps; r>0 && j < seq.length(); r--) {
					String readLetter = ""+seq.charAt(j++);
					String verLetter = ""+verSeq.charAt(verSeq.length()-r);
					
					if (!areTwoNucleotidesEqual(readLetter,verLetter)) {
						mm_encountered_here++;
						debugMes("walking read end gap: V[" + verLetter + "] vs. R[" + readLetter + "] ** conflict ** " , 20);
					}
					else {
						debugMes("walking read end gap: V[" + verLetter + "] vs. R[" + readLetter + "] OK " , 20);
					}
				}
				
			}
			else if (vertex_num_right_end_gaps > 0) {
				// over-ran the vertex sequence
				// shrink the sequence by the amount extended beyond the vertex
				j -= vertex_num_right_end_gaps;
				mm_encountered_here -= vertex_num_right_end_gaps;
				debugMes("because of vertex end gaps, walking read back by " + vertex_num_right_end_gaps + " bases.", 20);
				
				if (i != verSeq.length()) {
					debugMes("** ERROR: i=" + i + ", but verSeq.length() = " + verSeq.length(), 20);
				}
				
			}
			*/
			
			debugMes("mismatches encountered: " + mm_encountered_here, 20);
			
			
			if (mm_encountered_here >= zipper_mm && zipper_i == verSeq.length()) {
				debugMes("Zipper alignment mm: " + zipper_mm + " <= DP mm: " + mm_encountered_here +
						", so defaulting to earlier zipper alignment.", 20);
				i = zipper_i;
				j = zipper_j;
				mm_encountered_here = zipper_mm;
				max_left_gaps = 0;
				// retain any failed alignment status
			}
			else {
				
				failed_alignment = false; // reset it as needed given that DP was ok.
				
			}
			
			
			numMM = mm_encountered_here + totalNumMM;
			
		}
		
		
		// note, i and j are now 1 more than the array index due to the last i++,j++ of the loop above.


		float current_alignment_divergence = numMM / (float) j;
		debugMes("alignment divergence up to seq pos " + j + 
				" = mm: " + numMM +
				", div:" + current_alignment_divergence, 20);
		
		float local_vertex_alignment_divergence = mm_encountered_here / (float) i;
		debugMes("local vertex alignment divergence = " + mm_encountered_here + " / " + i + " = " + local_vertex_alignment_divergence, 20);
		
		
		// examine the alignment at this vertex to see if it passes our requirements
		if (i >= MIN_SEQ_LENGTH_TEST_DIVERGENCE && local_vertex_alignment_divergence >= MAX_READ_LOCAL_SEQ_DIVERGENCE) {
			failed_alignment = true;
			debugMes("local divergence exceeds max allowed: " + MAX_READ_LOCAL_SEQ_DIVERGENCE, 20);
		}
		
		if (max_left_gaps > MAX_LEFT_END_GAPS ) {
			failed_alignment = true;
		}
		
		if (
				
				// cumulative alignment stats up to and including vertex do not meet thresholds
				(current_alignment_divergence > MAX_READ_SEQ_DIVERGENCE || 
				numMM > MAX_MM_ALLOWED )
				
				|| 
				// alignment to current vertex fails
				failed_alignment
				)
		{
			debugMes("read "+readName+" has too many mismatches ("+numMM+") or too many left gaps (" + max_left_gaps + ")",19);

			if (failed_alignment) {
				// store it so we don't try again from this position in the sequence and at this vertex position.
				best_path_memoization.put(read_vertex_start_pos_token, null);
			}
			return(null); // go back and try alternative vertex if available

		} 
		else if (j==seq.length() || graph.getSuccessors(fromV) == null)
		{

			/////////////////////////////////////////
			// Reached end of the read being threaded (or ran out of vertices to explore)
			/////////////////////////////////////////

			if (graph.getSuccessors(fromV) == null) {
				// tack on unaligned terminus of sequence as mismatches
				mm_encountered_here += seq.length() - j;
			}
			
			// reached base case for recursion.
			
			debugMes("Reached end of read sequence.  Read" + readName + " with length: " + seq.length() 
					+ " and base [" + j + "] ends at position  [" + i + "] within node: " 
					+ fromV.getID() + " totaling " + mm_encountered_here + " mismatches. ", 19);

			Path_n_MM_count best_path = new Path_n_MM_count(fromV.getID(), mm_encountered_here, i, j);
			
			best_path_memoization.put(read_vertex_start_pos_token, best_path);
			
			return(new Path_n_MM_count(best_path));

		}
		else if (i==verSeq.length()) // move to the next ver
		{

			/////////////////////////
			// Reached end of vertex
			/////////////////////////


			// vertex sequence fully traversed, examine children vertices

			// Going on to recursive path mapping for read


			debugMes("Reached end of node sequence.  Read" + readName 
					+ " base [" + j + "] ends at position  [" + i + "] within node: " 
					+ fromV.getID() + " totaling " + mm_encountered_here + " mismatches. ", 19);


			// get list of the next vertices to explore.
			// examine edges prioritized by weight
		
		
			List<Integer> continueVersIds = new ArrayList<Integer>();
			List<SeqVertex> continueVers = new ArrayList<SeqVertex>();
			

			// just use simple successor vertices
			if (graph.getSuccessors(fromV) != null) {
				for (SeqVertex to_v : graph.getSuccessors(fromV)) {

					boolean check_initial_vertex_chars_match_seq = false;

					if (check_initial_vertex_chars_match_seq) {
						debugMes("-checking that next characters match up: " + to_v.getNameKmerAdj().charAt(0) + " vs. " + seq.charAt(j), 21);
						if (to_v.getNameKmerAdj().charAt(0)==seq.charAt(j)) {
							continueVers.add(to_v);
							continueVersIds.add(to_v.getID());


						}
					}
					else {
						continueVers.add(to_v);
						continueVersIds.add(to_v.getID());

					}

				}
			}
			
			debugMes("-reached end of vertex: " + fromV.getID() + ", exploring next vertices for continued path extension: " + continueVersIds, 19);

			
			Path_n_MM_count best_path = null;
			boolean tied_best = false;
			
			
			debugMes("Pursuing extension from : " + fromV.getShortSeqWID() + " to successors: " + continueVers, 19);
			
			List<Path_n_MM_count> all_best_paths_explored = new ArrayList<Path_n_MM_count>();
			
			for (Integer successor_vertex_id : continueVersIds) {
				
				
				debugMes("Exploring extension from node: " + fromV.getID() + " to node: " + successor_vertex_id, 20);
				
				Path_n_MM_count best_extension = updatePathRecursively(graph,successor_vertex_id,
																	seq, 
																	j,
																	KMER_SIZE-1, numMM, readName,
																	best_path_memoization);

				
				////////////////////////////////////////////
				// Back from Recursive Call.
				// Evaluate best paths from the successors.
				////////////////////////////////////////////

				
				/* testing for local sequence divergence within the alignment itself
				
				// first, check to see if it's an extension worth considering, given our local sequence divergence restrictions.
				if (best_extension != null
						&&
						best_extension.mismatch_count / (float) (seq.length() -1 - j) > MAX_READ_LOCAL_SEQ_DIVERGENCE) {
					
					debugMes("\tencountered max read sequence divergence: " + best_extension.mismatch_count / (float) (seq.length() -1 - j)
							+ ", disallowing extension: ." + best_extension, 19);
					
					best_extension = null;  // nullify the current best extension from successor_vertex_id
					
				}
				
				*/
				
				if (best_extension == null) {
					debugMes("\n\tFailed extension from " + fromV.getID() + " to : " + successor_vertex_id +  ".", 19);

				}
				else {
					// have a best extension

					all_best_paths_explored.add(best_extension);
					
					debugMes(readName + " best path so far from vertex: " + fromV.getID() 
							+ " to : " + successor_vertex_id
							+ " = " + best_extension.path +
							", with total mm: " + best_extension.mismatch_count, 20);
					
					if (best_path == null 
							||
							(best_extension.mismatch_count <= best_path.mismatch_count) ) {


						// test for tie condition
						if (best_path != null) {
							if (best_extension.mismatch_count == best_path.mismatch_count) {
								tied_best = true;
								debugMes("WARNING, Tied paths from vertex [V" + fromV_id + 
										" ]: \nPath A:\n" + best_extension + 
										"\nvs. Path B:\n" + best_path, 15);
							}
							else
								tied_best = false;
						}

						best_path = best_extension;

					}



				}
			} // end of successor vertex search.

			
			
			debugMes("Done with exploring paths from vertex: " + fromV.getID(), 20);
			debugMes("Paths and scores found are: ", 20);
			for (Path_n_MM_count pmm: all_best_paths_explored) {
				debugMes("\texplored path: " + pmm.path + " w/ mm: " + pmm.mismatch_count, 20);
			}
			if (best_path != null) {
				debugMes("\tAND best selected was: " + best_path.path + " w/ mm: " + best_path.mismatch_count, 20);
			}
			
			
			if (best_path != null) {
				
				
				if (tied_best) {
					debugMes("WARNING: TIED_READ_PATH", 15); 
					
					
					boolean TRUNCATE_TIED_PATH = false;
					
					if (TRUNCATE_TIED_PATH) {
						// truncate
						Path_n_MM_count truncated_path = new Path_n_MM_count(fromV.getID(), mm_encountered_here, i, j);

						best_path_memoization.put(read_vertex_start_pos_token, truncated_path);
						
						return(new Path_n_MM_count(truncated_path));
					}
					else {

						// add current node and local mismatches encountered here.
						best_path.add_path_n_mm(fromV.getID(), mm_encountered_here, i, j);

						best_path_memoization.put(read_vertex_start_pos_token, best_path);

						return(new Path_n_MM_count(best_path));

						
					}
				}
				else {
					// not a tie
					// add current node and local mismatches encountered here.
					
					best_path = new Path_n_MM_count(best_path);
					
					best_path.add_path_n_mm(fromV.getID(), mm_encountered_here, i, j);
					
					best_path_memoization.put(read_vertex_start_pos_token, best_path);
					
					return(new Path_n_MM_count(best_path));
				}
			}
			else {
				
				best_path_memoization.put(read_vertex_start_pos_token, null);
				return(null); // no extension possible.
			}
		}




		// should never end up here
		throw(new RuntimeException("should never end up here, supposedly.  i="+i 
					+ ", j=" + j + " ver length = "  + verSeq.length() + "  and readSeq length = " + seq.length() ));

	}


		
		
		
		

	/**
	 * create a hash that hold all the original vertices ids and the new ones 	
	 * @param graph
	 * @param rootIDs 
	 * @return the hash
	 */
	private static HashMap<Integer, LocInGraph> getOriginalVerIDsMappingHash(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {

		// clear double entries in the prevID list - not sure why they happen?
		for (SeqVertex v : graph.getVertices())
			v.clearDoubleEntriesToPrevIDs();


		HashMap<Integer, LocInGraph> hash = new HashMap<Integer,LocInGraph>();
		for (SeqVertex v : graph.getVertices())
		{

			// update the node tracker
			SeqVertex.nodeTracker.put(v.getID(), v); // beware - shouldn't have to do this, but finding myself having to due to getSeqVertex(id) not returning the correct vertex (with full sequence attached).
			
			debugMes("Graph vertex: " + v.getID() + " has sequence: " + v.getNameKmerAdj(), 22);
			
			Integer loc = 0;
			Integer vid = v.getID();
			// if the node id is new, than the real start is in the vector 
			if (vid>LAST_REAL_ID)
				loc = loc-1;
			else
			{
				debugMes("adding to "+vid+": Location of original node "+vid+" in index "+loc,20);
				hash.put(vid, new LocInGraph(vid,loc));

			}

			for (Vector<Integer> vec :  v.getPrevVerIDs())
			{
				loc++;
				for (Integer id : vec)
				{
					debugMes("adding to "+id+": Location of original node "+v.getID()+" in index "+loc,20);
					hash.put(id, new LocInGraph(v.getID(),loc));
				}
			}
		}
		return hash;
	}

	/**
	 * go over the graph file, and count the in flow and out flow of each node
	 * @param firstLetter 
	 * @throws IOException 
	 */
	private static void preProcessGraphFile(String filename,
			HashMap<Integer, Integer> outFlow, HashMap<Integer, Integer> inFlow, HashMap<Integer,String> kmers) throws IOException {

		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 
		String l = fileB.readLine(); // read header of component
		Integer from, to, supp;

		while (fileB.ready())

		{
			l = fileB.readLine();
			//	0       -1      3       ATTGAAAGCAAGTTTTCTTCGAAT        0
			//	1       0       3       TTGAAAGCAAGTTTTCTTCGAATT        0
			//	to		from	supp	kmer							stam       
			String[] fields = l.split("\t");
			from = Integer.parseInt(fields[1]);
			to = Integer.parseInt(fields[0]);
			supp = Integer.parseInt((fields[2]));
			String kmer = fields[3];

			if (!outFlow.containsKey(from))
				outFlow.put(from, supp);
			else
				outFlow.put(from, outFlow.get(from)+supp);

			if (!inFlow.containsKey(to))
				inFlow.put(to, supp);
			else
				inFlow.put(to, inFlow.get(to)+supp);

			kmers.put(to,kmer);
		}
	}


	/**
	 * given the filename, make a graph out of the connected components
	 * This time, keep the first letter of each kmer:
	 * keep the whole kmer, and then if there is an edge out, leave only first letter 
	 * @param filename
	 * @param rootIDs 
	 * @param inFlow in flow for all vertices
	 * @param outFlow out flow for all vertices
	 * @param firstLetter 
	 * @return
	 * @throws IOException
	 */
	private static DirectedSparseGraph<SeqVertex, SimpleEdge> buildNewGraphUseKmers(
			String filename, 
			Vector<Integer> rootIDs, HashMap<Integer,Integer> outFlow, 
			HashMap<Integer,Integer> inFlow, 
			HashMap<Integer, String> kmers) 
					throws IOException
					{

		BufferedReader fileB = 	new BufferedReader(new FileReader(filename)); 
		DirectedSparseGraph<SeqVertex, SimpleEdge> graph = 
				new DirectedSparseGraph<SeqVertex,SimpleEdge>();
		String l = fileB.readLine(); // read header of component
		Integer from, to;
		double supp;
		int linecount = 0;
		while (fileB.ready())

		{
			linecount++;
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 18 && linecount % 17 == 0) {
				System.err.print("\r[" + linecount + "]  ");
			}


			l = fileB.readLine();
			//	0       -1      3       ATTGAAAGCAAGTTTTCTTCGAAT        0
			//	1       0       3       TTGAAAGCAAGTTTTCTTCGAATT        0
			//	to		from	supp	kmer							stam       
			String[] fields = l.split("\t");
			from = Integer.parseInt(fields[1]);
			to = Integer.parseInt(fields[0]);
			supp = Double.parseDouble((fields[2]));
			if (supp < INITIAL_EDGE_ABS_THR )
				continue;

			
			// just tracking node ID values.
			if (from>LAST_ID)
				LAST_ID = from;

			if (to>LAST_ID)
				LAST_ID = to;

			String kmer = fields[3];
			if (KMER_SIZE == 0) {
				KMER_SIZE = kmer.length();
				debugMes("KMER_SIZE=" + KMER_SIZE, 5);
			}
			else if (KMER_SIZE != kmer.length()) {
				throw new RuntimeException("Error, discrepancy among kmer lengths.  Stored: " + KMER_SIZE + ", found: " + kmer.length() + "\n" + l);
				
			}
				
				
			
			SeqVertex fromV = getSeqVertex(graph, from);
			
			if (fromV==null && from>=0)
			{
				//fromV = new SeqVertex(from,firstLetter.get(from)+""+kmer.substring(0,KMER_SIZE-1));
				fromV = new SeqVertex(from,kmers.get(from));
				graph.addVertex(fromV);
			}
			
			boolean isRoot = (from<0 || fromV==null);

			
			SeqVertex toV = getSeqVertex(graph, to); // important to call this after possibly creating fromV, in case fromV == toV (otherwise, creating that vertex twice!) // bugfix
			
			if (isRoot)
			{
				if (toV==null)
				{
					toV = new SeqVertex(to, kmer, supp);
					graph.addVertex(toV);
					rootIDs.add(to);
				}
			}
			else
			{
				if (toV==null)
				{
					toV = new SeqVertex(to, kmer);
					graph.addVertex(toV);
				}
				SimpleEdge e = new SimpleEdge(supp, fromV.getID(), toV.getID()); 
				graph.addEdge(e, fromV, toV);
			}

		}

				
		return graph;
	}


	private static void writeDotFile(DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			String output_filename, String graphname, boolean printFullSeq) throws Exception
	{
		PrintStream p = new PrintStream(new FileOutputStream(output_filename));
		writeDotFile(graph,p,graphname,printFullSeq);
		p.close();
		
	}
	
	
	

	/**
	 * white to dot file with shortened seqs 
	 * @param graph
	 * @param p
	 * @param name
	 */
	private static void writeDotFile(DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			PrintStream p, String name)
	{
		writeDotFile(graph,p,name,false);
	}

	/**
	 * Write to dot file, where the list of paths are colored red -> blue
	 * @param graph
	 * @param p
	 * @param name
	 * @param vertices which vertices to print
	 */
	private static void writeDotFile(DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			PrintStream p, String name,boolean printFullSeq)
	{

		SeqVertex.set_graph(graph);
		
		//p.println("digraph "+name+"{");
		p.println("digraph G {");
		SeqVertex toVertex;
		int weight;

		//for each edge decide it's color
		for (SeqVertex vertex : graph.getVertices())
		{ //go over all vertices

			
			debugMes("Vertex: " + vertex.getShortSeqWconnectingIDs(graph), 15);
			
			
			String verDesc = ""+vertex.getID()+" [label=\"";
			if (printFullSeq)
				verDesc = verDesc.concat(""+vertex.getLongtSeqWID() + "[L:"+vertex.getNameKmerAdj().length()+"]\"");
			else {
				
				boolean include_discovery_time = true;
				if (include_discovery_time) {
					verDesc = verDesc.concat(""+vertex.getShortSeqWID() + "[L:"+vertex.getNameKmerAdj().length()+"]"
							+ "[T:" + vertex._dfsDiscoveryTime + "]\"");
				}
				else {
					verDesc = verDesc.concat(""+vertex.getShortSeqWID() + "[L:"+vertex.getNameKmerAdj().length()+"]\"");
				}
			}
			if (vertex.getWeightAvg()>25)
				verDesc = verDesc.concat(" ,style=bold,color=\"#AF0000\"");

			verDesc = verDesc.concat("]");
			if (!vertex.equals(T_VERTEX) && !vertex.equals(ROOT))
				p.println(verDesc);


			for (SimpleEdge edge : graph.getOutEdges(vertex)) //get all edges of vertex->?
			{
				toVertex = graph.getDest(edge);

				weight = (int) Math.round(edge.getWeight());
				String edgeStyle = "[label="+ weight +"]";

				if (weight>20)
					edgeStyle = "[style=bold,label="+ weight +",color=\"#AF0000\"]";
				if (!toVertex.equals(T_VERTEX) && !vertex.equals(ROOT))
					p.println(vertex.getID() + "->" + toVertex.getID() +edgeStyle);

			}
		}


		p.println("}");

	}



	/**
	 * Compact the given graph:
	 * for each vertex, if degree out = degree in = 1, and nextVertexIn ==1, remove this vertex, and connect edges
	 * @param graph
	 */
	private static boolean compactLinearPaths(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		debugMes("SECTION\n=================\nCOMPACTING THE GRAPH\n=================",5);
		//compact vertices
		Vector<SeqVertex> removeVertices = new Vector<SeqVertex>();
		Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();
		boolean changed = false;
		for (SeqVertex v1 : graph.getVertices())
		{
			//			debugMes("looking at vertex: "+v1);
			while (!v1.equals(ROOT) && graph.outDegree(v1)==1 )
			{
				// get out edge, should only have one.
				SimpleEdge e = null;
				for (SimpleEdge ei : graph.getOutEdges(v1))
					e = ei;
				
				// get the vertex attached.
				SeqVertex v2 = graph.getDest(e);
				if (graph.inDegree(v2)!=1 || v2.isToBeDeleted() || v2.equals(T_VERTEX) || v1.equals(v2)) {
					
					// avoiding loops and vertices that are in-branched.
					break;
				}
			
				
				
				debugMes("Found potential edge: "+e +" between "+v1 +" and "+v2,20);
				v1.concatVertex(v2, e.getWeight(),LAST_REAL_ID);
				debugMes("removing vertex "+v2+" was concatenated into "+v1,20);

				
				removeVertices.add(v2);
				v2.setToBeDeleted(true);

				changed = true;
				removeEdges.clear();
				for (SimpleEdge e2 : graph.getOutEdges(v2))
				{
					SeqVertex v3 = graph.getDest(e2);
					debugMes("Want to move edge " + e2 + "("+v2 +"->"+v3+") to ("+v1+"->"+v3,20);

					SimpleEdge newEdge = new SimpleEdge(e2, v1.getID(), v3.getID());
					graph.addEdge(newEdge, v1, v3);
					debugMes("\tadding edge: " + v1 + " to " + v3, 20);

					removeEdges.add(e2);
				}

				for (SimpleEdge re : removeEdges)
				{
					debugMes("removing edge " + re + "("+graph.getSource(re) +"->"+graph.getDest(re)+")",20);
					graph.removeEdge(re);
				}
				debugMes("removing edge " + e + "("+v1 +"->"+v2+")",20);
				graph.removeEdge(e);

			}
		}
		//remove all vertices that we don't want
		for (SeqVertex v : removeVertices)
		{
			graph.removeVertex(v);
		}
		return changed;
	}



	/**
	 * remove light edges from the graph. return true if something has changed
	 * @param graph
	 * @return
	 */
	private static boolean removeLightEdges(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		debugMes("removeLightEdges()", 10);

		boolean comp = false ; //removeLightCompEdges(graph);
		boolean in = removeLightInEdges(graph);
		boolean out = removeLightOutEdges(graph);
		boolean flow = removeLightFlowEdges(graph);
		return comp || in || out || flow;
	}


	/**
	 * Given a graph, go over all vertices and remove incoming or outgoing edges that do not match the flow (<2% coverage) see FLOW_THR
	 * When considering flow, this considers both the incoming and outgoing edges, but also the average node coverage.
	 * @param graph
	 * @return true if graph was changed.
	 */
	private static boolean removeLightFlowEdges(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		debugMes("SECTION\n=================\nREMOVING LIGHT FLOW EDGES\n=================",5);

		boolean changed = false;
		Collection<SeqVertex> all_vertices = graph.getVertices();
		
		int vertex_count = 0;
		int num_total_vertices = all_vertices.size();
		for (SeqVertex v : all_vertices)
		{
			
			vertex_count++;
			
			debugMes("Analyzing vertex: " + v.getID() + ", entry " + vertex_count + " of " + num_total_vertices, 25);
			int inDegree = graph.inDegree(v);
			int outDegree = graph.outDegree(v);
			debugMes("\thas inDegree: " + inDegree + ", outDegree: " + outDegree, 25);
			
		
			if (inDegree==0 && outDegree==0) {
				debugMes("\tSkipping vertex.", 25);
				continue;
			
			}
			
			int totalIn = 0, totalOut = 0;
			for (SimpleEdge e : graph.getInEdges(v))
				totalIn+=e.getWeight();

			for (SimpleEdge e : graph.getOutEdges(v))
				totalOut+=e.getWeight();

			debugMes("FLOW: total in for vertex "+v+" is "+totalIn + ", total out  is "+totalOut+", averageCov="+v.getWeightAvg(),20);

			Collection<SimpleEdge> removeEdges = new HashSet<SimpleEdge>();
			
			// out edges
			for (SimpleEdge e : graph.getOutEdges(v))
			{
				double e_avg_flow_thr_thresh = v.getWeightAvg() * FLOW_THR;
				if ( e.getWeight() < e_avg_flow_thr_thresh) {
					debugMes("EDGE_PRUNING::removeLightFlowEdges() removing low flow OUT edge " + e 
							+ " from "+ graph.getSource(e)+" to "+graph.getDest(e) +
							", FLOW_THR=" + FLOW_THR +
							", e_avg_flow_thr_thresh=: " + e_avg_flow_thr_thresh, 15);
					removeEdges.add(e);
				}
					
			}

			// in edges
			for (SimpleEdge e : graph.getInEdges(v))
			{
				double e_avg_flow_thr_thresh = v.getWeightAvg() * FLOW_THR;
				if (e.getWeight() < e_avg_flow_thr_thresh) {
					debugMes("EDGE_PRUNING::removeLightFlowEdges() removing low flow IN edge " + e 
							+ " from "+ graph.getSource(e)+" to "+graph.getDest(e) +
							", FLOW_THR=" + FLOW_THR +
							", e.weight=" + e.getWeight() + " < e_avg_flow_thr_thresh=: " + e_avg_flow_thr_thresh, 15);
					removeEdges.add(e);
				}
			}

			for (SimpleEdge re : removeEdges)
			{
				
				graph.removeEdge(re);
				changed = true;
			}
			
			debugMes("\tdone analyzing vertex: " + v.getID(), 25);
			
		}
		
		debugMes("== done removing Light Flow Edges.", 25);
		
		return changed;
	}
	/**
	 * go over the graph, and remove edges that are less than EDGE_THR (5%) from the rest of the entry flow
	 * @param graph
	 */
	private static boolean removeLightInEdges(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		debugMes("SECTION\n=================\nREMOVING LIGHT In EDGES\n=================",5);
		boolean somethingChanged = false;
		Queue<SeqVertex> allCurVers = new LinkedList<SeqVertex>(graph.getVertices());
		SeqVertex v = null;
		while ((v = allCurVers.poll())!=null)
		{
			if (graph.inDegree(v)<=1)
				continue;
			
			Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();
			
			// skip edges at simple cycles
			if (atSimpleCycle(graph, v)) {
				continue;
			}
			
			int totalIn = 0;
			for (SimpleEdge inE : graph.getInEdges(v))
			{
				totalIn+=inE.getWeight();
			}

			for (SimpleEdge inE : graph.getInEdges(v))
			{
				double e_edge_thr = totalIn*EDGE_THR;
				if (inE.getWeight() <= e_edge_thr)
				{
					debugMes("EDGE_PRUNING::removeLightInEdges() removing the edge: "+
							graph.getSource(inE)+" " + graph.getSource(inE).getName() + 
							" -> " +
							graph.getDest(inE)+ " " + graph.getDest(inE).getName() +
							" (weight: "+inE.getWeight()+" <= e_edge_thr: " + e_edge_thr +
							", EDGE_THR=" + EDGE_THR, 15); 
					removeEdges.add(inE);
					somethingChanged = true;
				}
			}
			for (SimpleEdge e : removeEdges)
				graph.removeEdge(e);
		}
		return somethingChanged;
	}


	private static boolean atSimpleCycle(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v) {
		for (SimpleEdge e : graph.getInEdges(v)) {
			if (graph.findEdge(v, graph.getSource(e)) != null)
				return(true);
		}
		for (SimpleEdge e: graph.getOutEdges(v)) {
			if (graph.findEdge(graph.getDest(e), v) != null)
				return(true);
		}
		
		return(false);
	}


	/**
	 * go over the graph, and remove edges that are less than EDGE_THR (10%) from the rest of the exit flow
	 * @param graph
	 */
	private static boolean removeLightOutEdges(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{
		debugMes("SECTION\n=================\nREMOVING LIGHT OUT EDGES\n=================",5);
		boolean somethingChanged = false;

		Queue<SeqVertex> allCurVers = new LinkedList<SeqVertex>(graph.getVertices());
		SeqVertex v = null;
		while ((v = allCurVers.poll())!=null)
		{
			if (graph.outDegree(v)<=1)
				continue;

			// skip edges at simple cycles
			if (atSimpleCycle(graph, v)) {
				continue;
			}
			
			Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();
			int totalOut = 0;
			for (SimpleEdge outE : graph.getOutEdges(v))
			{
				totalOut+=outE.getWeight();
			}

			for (SimpleEdge outE : graph.getOutEdges(v))
			{
				double e_edge_thr = totalOut * EDGE_THR;
				if (outE.getWeight() <= e_edge_thr)
				{
					
					debugMes("EDGE_PRUNING::removeLightOutEdges() removing the edge: " +
							graph.getSource(outE)+ " " + graph.getSource(outE).getName() +
							" -> " + 
							graph.getDest(outE)+ " " + graph.getDest(outE).getName() +
							" (weight: "+outE.getWeight()+" <= e_edge_thr: " + e_edge_thr +
							", EDGE_THR=" + EDGE_THR, 15);
					
					removeEdges.add(outE);
					somethingChanged = true;
				}
			}
			for (SimpleEdge e : removeEdges)
				graph.removeEdge(e);

		}
		return somethingChanged;
	}

	/**
	 * Return the SeqVertex with the given id within the given graph.
	 * @param graph
	 * @param id
	 * @return
	 */
	private static SeqVertex getSeqVertex(DirectedSparseGraph<SeqVertex, SimpleEdge> graph, int id)
	{

		SeqVertex v = SeqVertex.retrieveSeqVertexByID(id);
		
		if (graph.containsVertex(v)) {
			return(v);
		}
		else {
			return(null);  // note, SeqVertex stores all vertices even if removed from the graph, so let's ensure it's still in the graph.
		}

		/*  orig code too slow for large graphs
		for (SeqVertex v : graph.getVertices())
		{
			if (v.getID() == id)
				return v;
		}
		return null;
		 */

	}


	/**
	 * Given the string seq, return it in fasta format
	 * @param seq - seq
	 * @param name - seq name
	 * @return
	 */
	private static String getSeqFasta(String seq,String name){
		String res = "";
		res = res.concat(">"+name+"\n");

		int i=0;
		for (; i<seq.length()-LINE_LEN ; i+=LINE_LEN)
		{
			res = res.concat(seq.substring(i, i+LINE_LEN)+"\n");
		}
		res = res.concat(seq.substring(i)+"\n");
		return res;
	}


	/**
	 * return the next available vertex id.
	 * @return
	 */
	private static int getNextID() {
		LAST_ID++;
		return LAST_ID;
	}

	/**
	 * return a topological order on the graph's vertices.
	 * @param graph
	 * @return list of nodes.
	 */
	private static List<SeqVertex> getTopologicalOrder(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, My_DFS dfs) {

		
		Map<SeqVertex,Number> finished = dfs.getFinishing();
		SeqVertexFinishTimeComparator finishingTimeComparator = new SeqVertexFinishTimeComparator();
		debugMes("getTopologicalOrder(), Vertex count: " + graph.getVertexCount(), 18);


		PriorityQueue<SeqVertex> fQueue = new PriorityQueue<SeqVertex>(graph.getVertexCount(),finishingTimeComparator  );

		for (SeqVertex v : finished.keySet())
		{
			fQueue.add(v);
		}
		List<SeqVertex> order = new ArrayList<SeqVertex>();
		while (!fQueue.isEmpty())
		{
			order.add(fQueue.poll());
		}
		return order;
	}



	/**
	 * Go over each sub component of the given graph, and calc the following:
	 * total coverage (sum of weights)
	 * average coverage
	 * number of paths
	 * @param graph
	 */
	private static void calcSubComponentsStats(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {

		Set<Set<SeqVertex>> comps = divideIntoComponents(graph);
		int numComp = comps.size();
		for (Set<SeqVertex> comp : comps)
		{
			//now we have one comp in hand
			Vector<Double> allW = new Vector<Double>();
			int compID = -1;
			for (SeqVertex v : comp)
			{
				if (compID==-1)
					compID = v.getID();

				allW.addAll(0, v.getWeights());
				for (SimpleEdge outE : graph.getOutEdges(v))
				{
					allW.add(0, outE.getWeight());
				}
			}
			SeqVertex v1 = getSeqVertex(graph, compID);
			if (allW.size()==0 || (comp.size()==1 && v1.getName().length()<MIN_OUTPUT_SEQ))
			{
				//this is a single node with a single letter
				debugMes("removing component with node "+compID+" which has only one node with short seq "+v1.getName(),20);
				graph.removeVertex(v1);
				numComp = numComp-1;
				continue;

			}

			int t=0;
			for (Double w: allW)
				t+=w;

			//System.err.println("t=" + t + ", allW.size()=" + allW.size());
			float avgCov = (float)t/allW.size();
			debugMes("SubComp: "+compID+" has "+ comp.size() +" nodes; total coverage: "+t+" average: "+avgCov,20);

			if (avgCov<COMP_AVG_COV_THR-0.5) //FIXME: added 0.5 for testing with low cov seq
			{
				debugMes("removing component with node "+compID+" which has only average coverage of "+
						avgCov+ " < "+COMP_AVG_COV_THR,20);
				for (SeqVertex v : comp)
					graph.removeVertex(v);
				numComp = numComp-1;
			}
		}
		debugMes("number of good components: "+numComp,10);
	}

	/**
	 * divide the graph into its components
	 * @param graph
	 * @return set of components
	 */
	private static Set<Set<SeqVertex>> divideIntoComponents(DirectedSparseGraph<SeqVertex, SimpleEdge> graph) 
	{

		WeakComponentClusterer<SeqVertex, SimpleEdge> compClus = new WeakComponentClusterer<SeqVertex, SimpleEdge>();
		Set<Set<SeqVertex>> comps = compClus.transform(graph);
		
		int comp_counter = 0;
		for (Set<SeqVertex> s : comps) {
			debugMes("\nComponentDivision: " + comp_counter + " contains the following vertices:", 10);
			for (SeqVertex v : s) {
				debugMes("node_id: " + v.getID(), 10);
			}
			comp_counter++;
			
		}
		
		return comps;

	}

	/**
	 * connect the source node to each node with indegree=0,
	 * connect each node with outdegree=0 to the target node 
	 * Also add reads from the root to each of the nodes, and from the ends too.
	 * @param graph
	 * @param comp the current component
	 * @param combinedReadHash 
	 */
	private static void addSandT(DirectedSparseGraph<SeqVertex, SimpleEdge> graph, Set<SeqVertex> comp, HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash)
	{
		//		debugMes("=================\nADDING S AND T\n=================",10);

		graph.addVertex(ROOT);
		graph.addVertex(T_VERTEX);
		SimpleEdge e=null;
		//		for (SeqVertex v : graph.getVertices())
		for (SeqVertex v : comp)
		{
			if (graph.inDegree(v)==0 && !v.equals(ROOT) && !v.equals(T_VERTEX)) // connect S to this vertex
			{
				double w = v.getFirstWeight();
				if (w==-1) // single letter node?
				{
					debugMes("got a single letter node here.. "+v,20);
					w = 1;
				}
				e = new SimpleEdge(w, ROOT.getID(), v.getID());
				graph.addEdge(e, ROOT, v);

				debugMes("Adding edge from S to "+v,20);

				PairPath pathD = new PairPath();
				pathD.addToPath1(ROOT.getID());
				pathD.addToPath1(v.getID());
				
				if (!combinedReadHash.containsKey(ROOT.getID()))
					combinedReadHash.put(ROOT.getID(), new HashMap<PairPath,Integer>());
				combinedReadHash.get(ROOT.getID()).put(pathD, MIN_READ_SUPPORT_THR);
				
				
				/*
				for (SeqVertex v2 : graph.getSuccessors(v))
				{
					PairPath pathD = new PairPath();
					pathD.addToPath1(ROOT.getID());
					pathD.addToPath1(v.getID());
					pathD.addToPath1(v2.getID());
					if (!combinedReadHash.containsKey(ROOT.getID()))
						combinedReadHash.put(ROOT.getID(), new HashMap<PairPath,Integer>());
					combinedReadHash.get(ROOT.getID()).put(pathD, MIN_READ_SUPPORT_THR);

				}		
				*/		
			}

			if (graph.outDegree(v)==0 && !v.equals(T_VERTEX) && !v.equals(ROOT)) // connect this vertex to T
			{
				double w = v.getLastWeight();
				if (w==-1)
					w=1;
				e = new SimpleEdge(w, v.getID(), T_VERTEX.getID());
				graph.addEdge(e, v, T_VERTEX);
				debugMes("Adding edge from "+v+" to T",20);


				PairPath pathD = new PairPath();
				pathD.addToPath1(v.getID());
				pathD.addToPath1(T_VERTEX.getID());
				
				if (!combinedReadHash.containsKey(v.getID()))
					combinedReadHash.put(v.getID(), new HashMap<PairPath,Integer>());
				combinedReadHash.get(v.getID()).put(pathD, MIN_READ_SUPPORT_THR);
				
				/*
				for (SeqVertex v2 : graph.getPredecessors(v))
				{
					PairPath pathD = new PairPath();
					pathD.addToPath1(v2.getID());
					pathD.addToPath1(v.getID());
					pathD.addToPath1(T_VERTEX.getID());
					if (!combinedReadHash.containsKey(v2.getID()))
						combinedReadHash.put(v2.getID(), new HashMap<PairPath,Integer>());
					combinedReadHash.get(v2.getID()).put(pathD, MIN_READ_SUPPORT_THR);


				}
				*/
			}
		}
	}

	/**
	 * given the graph, remove all edges of S and T
	 * @param graph
	 */
	private static void removeAllEdgesOfSandT(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		Set<SimpleEdge> removeEdges = new HashSet<SimpleEdge>();
		if (graph.containsVertex(ROOT))
			for (SimpleEdge e : graph.getOutEdges(ROOT))
				removeEdges.add(e);
		if (graph.containsVertex(T_VERTEX))
			for (SimpleEdge e : graph.getInEdges(T_VERTEX))
				removeEdges.add(e);

		for (SimpleEdge re : removeEdges)
			graph.removeEdge(re);

	}

	 
	private static boolean dealWithLoops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
			Set<SeqVertex> comp, 
			HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash) {

		boolean res = false;
		DijkstraShortestPath<SeqVertex, SimpleEdge> dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(graph);
		
		///////////////////////////////////////////////////////////////////////////////////////
		//  ===== Capture all remaining loops found. ========== 
		//
		// These should be only those repeats that aren't evident in the individual read paths,
		// since the read-evident repeats were unrolled earlier.
		/////////////////////////////////////////////////////////////////////////////////////////
		
		Set<Set<SimpleEdge>> curLoops = new HashSet<Set<SimpleEdge>>();

		// find all loops in the graph by seeing if, given edge v->v2, there is a path from v2 back to v
		for (SeqVertex v : comp)
		{
			for (SeqVertex v2 : graph.getSuccessors(v))
			{
				if (dp.getDistance(v2, v)!=null) // there is a connection between v->v2->... ->v
				{
					//path has all edges from v to itself thru v2
					List<SimpleEdge> loopPath = dp.getPath(v2, v);
					
					// v2 is successor of v, so let's just add the v->v2 edge too, complete the full loop.
					loopPath.add(0, graph.findEdge(v, v2)); 
					
					// capture the path IDs for debugMes reporting below.
					List<Integer> pathIDs = new ArrayList<Integer>();
					for (SimpleEdge e : loopPath)
						pathIDs.add(graph.getDest(e).getID());

					// Collect the loop edge set.
					Set<SimpleEdge> loopPath_set = new HashSet<SimpleEdge>(loopPath);
				
					if (!curLoops.contains(loopPath_set))
					{
						curLoops.add(loopPath_set);
						debugMes("adding the loop path "+pathIDs+" to the curLoops",12);
					}else
					{
						debugMes("not adding the loop path "+pathIDs+" to the curLoops",12);
					}
				}
			}
		}

		if (curLoops.isEmpty())
			return res;


		////////////////////////
		// process found loops
		///////////////////////


		Set<SimpleEdge> allRelevantEdges = new HashSet<SimpleEdge>();

		for (Set<SimpleEdge> loopPath_set : curLoops)
			for (SimpleEdge e : loopPath_set)
			{
				e.increaseNumOfLoopsBy1();
				allRelevantEdges.add(e);
			}


		////////////////////////////////////////
		//  break complex loops
		////////////////////////////////////////

		
		if (!allRelevantEdges.isEmpty()){
			
			Comparator<Object> numLoopsComparator = new numLoopsEdgeComparator(graph);
			
			PriorityQueue<SimpleEdge> edgesQ = new PriorityQueue<SimpleEdge>(allRelevantEdges.size(), numLoopsComparator);
			
			edgesQ.addAll(allRelevantEdges);

			//while there are still loops
			// find the next edge that can be removed to reduce the number of loops
			// updated queue: remove all edges, and update their loop content
			
			SimpleEdge nextEtoRemove;
			while ( (!curLoops.isEmpty()) && (! edgesQ.isEmpty()) ) 
			{

				//FIXME: there was a situation where curLoops was not empty, 
				//       but edgesQ was, so I added edgesQ to the while condition. Investigate why this might happen.
				//       In this case, a node was involved in a self loop and a double-loop.

				nextEtoRemove = edgesQ.poll();

				if (graph.getSource(nextEtoRemove) == null 
						|| graph.getDest(nextEtoRemove) == null 
						|| nextEtoRemove.getNumOfLoopsInvolved() <= 0) { 
					continue; 
				}

				debugMes("removing the edge " + graph.getSource(nextEtoRemove).getID() + "->" + 
						graph.getDest(nextEtoRemove).getID() + " that appears in "
						+nextEtoRemove.getNumOfLoopsInvolved() + " loops",15);
				
				
				// remove the loops that have this edge from curLoops
				Set<Set<SimpleEdge>> removeLoops = new HashSet<Set<SimpleEdge>>();
				for (Set<SimpleEdge> loopPath_set : curLoops)
					if (loopPath_set.contains(nextEtoRemove))
					{
						debugMes("the loop "+ loopPath_set+" is now solved",15);
						removeLoops.add(loopPath_set);

						// update the number of loops involved in each edge
						for (SimpleEdge e : loopPath_set)
							e.decreaseNumOfLoopsBy1();
					}
				
				for (Set<SimpleEdge> loopPath_set : removeLoops)
					curLoops.remove(loopPath_set);

				
				//update the queue. remove all, and insert again if numLoops>0.
				SimpleEdge[] relEdges = (SimpleEdge[]) edgesQ.toArray(new SimpleEdge[0]);
				edgesQ.clear();
				for (SimpleEdge otherE : relEdges)
					if (otherE.getNumOfLoopsInvolved()>0)
						edgesQ.add(otherE);

				// remove this edge
				graph.removeEdge(nextEtoRemove);
				res = true;
			}
		}
		return res;


	}


	
	
	
	/**
	 * given the graph and the node with the self loop,
	 * find the reads that support this loop, and multiply this vertex as many times as needed, and then remap these reads.
	 * @param graph
	 * @param v
	 * @param combinedReadHash
	 * @param newVers 
	 */
	private static void dealWithSelfLoops(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v,
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash, Set<SeqVertex> newVers) {


		int vid = v.getID();
		int maxNumOfOccurrences = 0;

		HashMap<PairPath, Integer> relaventReads = new HashMap<PairPath, Integer>();
		for (Integer startV : combinedReadHash.keySet())
		{
			for (PairPath path: combinedReadHash.get(startV).keySet())
			{
				int numOcc = path.numOccurrences(vid);
				if (numOcc>0)
				{					
					Integer count = combinedReadHash.get(startV).get(path);
					if (count == null)
						debugMes("stop here",10);
					relaventReads.put(path,count);
				}

				if ( numOcc> maxNumOfOccurrences) //this read includes this vertex
				{
					debugMes("the read "+path+" includes the vertex "+vid+" "+numOcc+" times",19);
					maxNumOfOccurrences = numOcc;
				}
			}
		}

		// remove the self loop
		SimpleEdge removeE = graph.findEdge(v, v);
		double oldW = removeE.getWeight();
		List<Integer> newVerIDs = new ArrayList<Integer>();
		newVerIDs.add(vid);

		graph.removeEdge(removeE);
		debugMes("removing the edge between "+ v +" and itself",20);
		// multiply this node maxNumOfOccurrences times
		int upID = vid;
		int downID = -1;
		ArrayList<SimpleEdge> removeEdges = new ArrayList<SimpleEdge>();
		for (int i=2; i<=maxNumOfOccurrences; i++)
		{
			if (downID!=-1)
				upID = downID;

			downID = getNextID();
			newVerIDs.add(downID);

			SeqVertex newV = new SeqVertex(downID, v);
			debugMes("adding the new vertex "+newV.getID(),20);

			newV.setOrigButterflyID(v.getID());
			graph.addVertex(newV);
			SeqVertex upV = getSeqVertex(graph, upID);

			newVers.add(newV);

			for (SeqVertex vOut : graph.getSuccessors(upV))
			{
				debugMes("adding an edge between "+newV.getID()+" and "+vOut.getID(),20);
				graph.addEdge(new SimpleEdge(graph.findEdge(v, vOut), newV.getID(), vOut.getID()), newV, vOut);
			}

			debugMes("adding an edge between "+upID+" and "+newV.getID(),20);
			graph.addEdge(new SimpleEdge(oldW, upV.getID(), newV.getID()), upV, newV);

		}
		// moved to the end of loop added in the new loop opening process - Feb 2013
		SeqVertex newV = getSeqVertex(graph, downID);
		for (SeqVertex vOut : graph.getSuccessors(v))
		{
			if (!newVerIDs.contains(vOut.getID())){
				debugMes("adding an edge between "+newV.getID()+" and "+vOut.getID(),0);
				SimpleEdge e = graph.findEdge(v, vOut);
				graph.addEdge(new SimpleEdge(e, newV.getID(), vOut.getID()), newV, vOut);
				debugMes("removing the edge between "+ v.getID() +" and "+vOut.getID(),20);
				removeEdges.add(e);
			}
		}

		//remove edges:
		for (SimpleEdge re : removeEdges){
			graph.removeEdge(re);

		}

		List<Integer> loopVIDs = new ArrayList<Integer>();
		loopVIDs.add(vid);
		List<List<Integer>> newVerIDsList = new ArrayList<List<Integer>>();
		newVerIDsList.add(newVerIDs);
		updateReadsAfterLoopOpening(combinedReadHash,relaventReads,loopVIDs,newVerIDsList,maxNumOfOccurrences);

	}

	/**
	 * Given the combinedReadHash, and the relevant reads, update their paths.
	 * @param combinedReadHash
	 * @param relevantReads
	 * @param loopVIDs
	 * @param newVerIDs
	 * @param maxNumOfOccurrences
	 */
	private static void updateReadsAfterLoopOpening(
			HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash,
			HashMap<PairPath, Integer> relevantReads, List<Integer> loopVIDs,
			List<List<Integer>> newVerIDs, int maxNumOfOccurrences) {

		for (PairPath path: relevantReads.keySet())
		{
			Integer origFirstV = path.getFirstID();
			Integer origCount = combinedReadHash.get(origFirstV).get(path);
			List<Integer> newPath1 = new ArrayList<Integer>(path.getPath1());
			List<Integer> newPath2 = new ArrayList<Integer>(path.getPath2());

			if (loopVIDs.size()==1)
			{
				updatePathOfSelfLoop(newPath1,loopVIDs,newVerIDs.get(0),maxNumOfOccurrences);
				updatePathOfSelfLoop(newPath2,loopVIDs,newVerIDs.get(0),maxNumOfOccurrences);
			} else
			{
				updatePathOfDoubleLoop(newPath1,loopVIDs,newVerIDs.get(0),newVerIDs.get(1),maxNumOfOccurrences);
				updatePathOfDoubleLoop(newPath2,loopVIDs,newVerIDs.get(0),newVerIDs.get(1),maxNumOfOccurrences);
			}
			// path hasn't changed
			if (path.getPath1().equals(newPath1) && path.getPath2().equals(newPath2))
				continue;

			// both are empty now
			if (newPath1.isEmpty() && newPath2.isEmpty())
				combinedReadHash.get(origFirstV).remove(path);

			// at least one has changed
			PairPath newKey;
			if (newPath1.isEmpty())	
				newKey = new PairPath(newPath2,new ArrayList<Integer>());
			else if (newPath2.isEmpty())
				newKey = new PairPath(newPath1,new ArrayList<Integer>());
			else
				newKey = new PairPath(newPath1,newPath2);


			Integer firstV = newKey.getFirstID();
			if (!combinedReadHash.containsKey(firstV))
				combinedReadHash.put(firstV, new HashMap<PairPath, Integer>());

			if (combinedReadHash.get(firstV).containsKey(newKey))
			{
				Integer oldCount = combinedReadHash.get(firstV).get(newKey);
				combinedReadHash.get(firstV).put(newKey,oldCount+origCount);
				combinedReadHash.get(firstV).remove(path);
			}else
			{
				combinedReadHash.get(firstV).put(newKey,origCount);
			}
			
			// remove the old loop-containing path
			combinedReadHash.get(origFirstV).remove(path);
		}

	}



	/**
	 * given a path, the vid of the self loop, and the new vertices' id, update the path
	 * if the path starts of ends inside the loop, trim this part of the path, and leave only the outside info.
	 * @param path
	 * @param vid
	 * @param newVerIDs
	 * @return
	 */
	private static void updatePathOfSelfLoop(List<Integer> path, List<Integer> loopVIDs,
			List<Integer> newVerIDs,int maxNumOcc) {
		int vid = loopVIDs.get(0).intValue();
		String origPath = ""+path;
		Set<Integer> loopVs = new HashSet<Integer>();
		loopVs.add(vid);
		boolean changed = false;
		if (path.contains(vid))
		{
			if (path.get(0).intValue()==vid)
			{ //starts inside the loop
				changed = true;
				if (path.get(path.size()-1).intValue()==vid)
					//starts and ends inside the loop
					if (path.size()==maxNumOcc)
					{
						for (int i=1 ; i<=path.size()-1 ; i++)
							path.set(i,newVerIDs.get(i));

						changed = true;
					}else
						path.clear();
				else
					updatePathToRemoveLoopNodes(path,loopVs);
			}else
			{ // starts and ends outside the loop
				for (int i=1 ; i<=path.size()-1 ; i++)
				{
					if (path.get(i).intValue()==vid) // i>0
					{
						int j = newVerIDs.indexOf(path.get(i-1));
						if (j>=0)
						{
							path.set(i, newVerIDs.get(j+1));
							changed = true;
						}
					}

				}
			}
		}
		if (changed)
			debugMes("path changed from "+origPath+" to "+path,20);
	}

	/**
	 * remove the integers that are inside the loop
	 * @param path
	 * @param loopVs
	 */
	private static void updatePathToRemoveLoopNodes(List<Integer> path,
			Set<Integer> loopVs) {
		List<Integer> indicesToRemove = new ArrayList<Integer>();
		for (int i=0 ; i<=path.size()-1 ; i++)
			if (loopVs.contains(path.get(i)))
				indicesToRemove.add(i);
		Collections.sort(indicesToRemove);
		Collections.reverse(indicesToRemove);
		for (Integer i : indicesToRemove)
			path.remove(i.intValue());
	}


	 /**
     * given the graph and the node with the self loop,
     * find the reads that support this loop, and multiply this vertex as many times as needed
     * @param graph
     * @param t_v1
     * @param t_v2
     * @param combinedReadHash
     * @param newVers 
     */
    private static void dealWithDoubleLoops(
                    DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex t_v1,
                    SeqVertex t_v2,
                    HashMap<Integer, HashMap<PairPath, Integer>> combinedReadHash, Set<SeqVertex> newVers) 
    {
            int v1_id=-1; // the one inside the regular flow
            int v2_id=-1; // the addition
            
            //CONTINUE HERE!

            if (graph.getSuccessorCount(t_v1)>1)
            {
                    v1_id = t_v1.getID();
                    v2_id = t_v2.getID();
            } else if (graph.getSuccessorCount(t_v2)>1)
            {
                    v1_id = t_v2.getID();
                    v2_id = t_v1.getID();
            }

            if (v1_id==-1)
            {               
                    //FIXME - decide randomly, doesn't solve the loops right. (ignores input edges to t_v1)
                    v1_id = t_v1.getID();
                    v2_id = t_v2.getID();
            }

            debugMes("\n\ndealWithDoubleLoops: v1 = " + v1_id + ", v2 = " + v2_id, 15);

            //////////////////////////////////////////////////////
            // count the number of times v2 appears within reads
            
            HashMap<PairPath, Integer> relevantReads = new HashMap<PairPath, Integer>();

            int maxNumOfOccurrences = 0;
            for (Integer startV : combinedReadHash.keySet())
            {
                    for (PairPath path: combinedReadHash.get(startV).keySet())
                    {
                            int numOcc2 = path.numOccurrences(v2_id);
                            debugMes("read-to-vertex-count: " + path + " contains vertex: " + v2_id + " this many times: " + numOcc2, 15);
                            if (numOcc2>0)
                            {
                                Integer count = combinedReadHash.get(startV).get(path);
                                if (count == null)
                                {
                                	// shouldn't happen

                                        for (PairPath path2: combinedReadHash.get(startV).keySet())
                                        {
                                                debugMes("path: "+path2+" with hashCode "+path2.hashCode(),15);
                                                debugMes("path: "+path2+" with value "+combinedReadHash.get(startV).get(path2),15);
                                        }
                                }
                                relevantReads.put(path,count);
                                
                                if ( numOcc2> maxNumOfOccurrences) //this read includes this vertex
                                {
                                        debugMes("MAX:the read "+path+" includes the vertex "+v2_id+" "+numOcc2+" times",15);
                                        maxNumOfOccurrences = numOcc2;
                                }

                        }
                }
        }
            
        
            
        // the loop is v1 (v2,v1)*
        //if we count how many v2 appears, then the number appearances of v1 is one more.
        SeqVertex v1 = getSeqVertex(graph, v1_id);
        SeqVertex v2 = getSeqVertex(graph, v2_id);
        
        
        List<Integer> newVerIDs_v1 = new ArrayList<Integer>();
        List<Integer> newVerIDs_v2 = new ArrayList<Integer>();
        newVerIDs_v1.add(v1_id);
        newVerIDs_v2.add(v2_id);

        // remove the self loop
        SimpleEdge removeE = graph.findEdge(v2, v1);
        double oldW = removeE.getWeight();
        double oldW2 = graph.findEdge(v1, v2).getWeight();

        graph.removeEdge(removeE);
        debugMes("removing the edge between "+ v2_id +" and "+v1_id,15);
        // multiply this node maxNumOfOccurrences times

        int up_v1 = v1_id;
        if (maxNumOfOccurrences>=1) //multiply only v1
        {
                SeqVertex newV = new SeqVertex(getNextID(), v1);
                newV.setOrigButterflyID(v1_id);
                graph.addVertex(newV);
                newVerIDs_v1.add(newV.getID());

                newVers.add(newV);

                //                      // removed in the new loop opening process - Feb 2013
                //                      for (SeqVertex vOut : graph.getSuccessors(v1))
                //                      {
                //                              if (!vOut.equals(v2))
                //                              {
                //                                      debugMes("adding an edge between "+newV.getID()+" and "+vOut.getID(),20);
                //                                      graph.addEdge(new SimpleEdge(graph.findEdge(v1, vOut)), newV, vOut);
                //                              }
                //                      }
                debugMes("adding an edge between "+v2_id+" and "+newV.getID(),20); 
                graph.addEdge(new SimpleEdge(oldW, v2.getID(), newV.getID()), v2, newV);

                up_v1 = newV.getID();
        }
        
        int up_v2 = v2_id;
        //int down_v1 = -1;
        int down_v2 = -1;

        int down_v1 = up_v1;
        
       
        for (int i=2; i<=maxNumOfOccurrences; i++) // multiple v2-v1
        {
                if (down_v1!=-1)
                {
                	// on next iteration
                    up_v1 = down_v1;
                }

                down_v1 = getNextID();
                down_v2 = getNextID();
                newVerIDs_v1.add(down_v1);
                newVerIDs_v2.add(down_v2);

                SeqVertex newV1 = new SeqVertex(down_v1, v1);
                newV1.setOrigButterflyID(v1_id);
                SeqVertex newV2 = new SeqVertex(down_v2, v2);
                newV2.setOrigButterflyID(v2_id);
                //debugMes("i="+i+"("+maxNumOfOccurrences+") adding newV1:"+newV1+" newV2:"+newV2,10);
                graph.addVertex(newV1);
                graph.addVertex(newV2);
                newVers.add(newV1);
                newVers.add(newV2);

                SeqVertex upV = getSeqVertex(graph, up_v1);
                
                /*
                SeqVertex orig_upV = getSeqVertex(graph, v1_id);

                if (i==maxNumOfOccurrences){// this cirteria added in the new loop opening process - Feb 2013  // moved down below as bugfix, since need to do this for no-repeats too.  Sept 2014
                        for (SeqVertex vOut : graph.getSuccessors(orig_upV))
                        {
                                if (!newVerIDs_v2.contains(vOut.getID()))
                                {
                                        debugMes("adding an edge between "+newV1.getID()+" and "+vOut.getID(),20);
                                        SimpleEdge e = graph.findEdge(orig_upV, vOut);
                                        graph.addEdge(new SimpleEdge(e), newV1, vOut);
                                        debugMes("removing an edge between "+orig_upV.getID()+" and "+vOut.getID(),20); 
                                        removeEdges.add(e);
                                }
                        }
                }
                
                */
                
                
                //                      // removed in the new loop opening process - Feb 2013
                //                      for (SeqVertex vIn : graph.getPredecessors(getSeqVertex(graph, up_v2)))
                //                      {
                //                              if (!newVerIDs_v1.contains(vIn.getID()))
                //                              {
                //                                      debugMes("$$adding an edge between "+vIn.getID()+" and "+down_v2,10);
                //                                      graph.addEdge(new SimpleEdge(graph.findEdge(vIn, getSeqVertex(graph, up_v2))), vIn, newV2);
                //                              }
                //                      }
                debugMes("adding an edge between "+up_v1+" and "+newV2.getID(),20);
                graph.addEdge(new SimpleEdge(oldW, upV.getID(), newV2.getID()), upV, newV2);

                debugMes("adding an edge between "+newV2.getID()+" and "+newV1.getID(),20);
                graph.addEdge(new SimpleEdge(oldW2, newV2.getID(), newV1.getID()), newV2, newV1);


        }

        List<Integer> loopVIDs = new ArrayList<Integer>();
        loopVIDs.add(v1_id);
        loopVIDs.add(v2_id);
        List<List<Integer>> newVerIDs = new ArrayList<List<Integer>>();
        newVerIDs.add(newVerIDs_v1);
        newVerIDs.add(newVerIDs_v2);

        // relocate original V1's out-edges
        SeqVertex orig_upV = getSeqVertex(graph, v1_id);
        SeqVertex down_v1_vertex = getSeqVertex(graph, down_v1);
   
        ArrayList<SimpleEdge> removeEdges = new ArrayList<SimpleEdge>();
        
        for (SeqVertex vOut : graph.getSuccessors(orig_upV))
        {
        	if (!newVerIDs_v2.contains(vOut.getID()))
        	{
        		debugMes("adding an edge between "+down_v1_vertex.getID()+" and "+vOut.getID(),20);
        		SimpleEdge e = graph.findEdge(orig_upV, vOut);
        		graph.addEdge(new SimpleEdge(e, down_v1_vertex.getID(), vOut.getID()), down_v1_vertex, vOut);
        		debugMes("removing an edge between "+orig_upV.getID()+" and "+vOut.getID(),20); 
        		removeEdges.add(e);
        	}
        }
        
        
        //remove edges:
        for (SimpleEdge re : removeEdges){
                graph.removeEdge(re);

        }
        updateReadsAfterLoopOpening(combinedReadHash,relevantReads,loopVIDs,newVerIDs,maxNumOfOccurrences);


    }



	/**
	 * given a path, the vid of the loop vertices, and the new vertices' id, update the path
	 * if the path starts of ends inside the loop, trim this part of the path, and leave only the outside info.
	 * @param path
	 * @param loopVIDs
	 * @param newVerIDsV1
	 * @param newVerIDsV2
	 * @param maxNumOfOccurrences
	 */
	private static void updatePathOfDoubleLoop(
			List<Integer> path, 
			List<Integer> loopVIDs,
			List<Integer> newVerIDsV1, 
			List<Integer> newVerIDsV2, 
			int maxNumOfOccurrences) {
		
		int v1_id = loopVIDs.get(0).intValue();
		int v2_id = loopVIDs.get(1).intValue();

		
		debugMes("updatePathOfDoubleLoop, read: " + path + ", loop vertices: " + v1_id + ", " + v2_id, 15);
		
		
		if (path.isEmpty())
			return;
		boolean changed = false;

		String origPath = ""+path;

		Set<Integer> loopVs = new HashSet<Integer>();
		loopVs.add(v1_id);
		loopVs.add(v2_id);

		int firstV = path.get(0).intValue();
		int lastV = path.get(path.size()-1).intValue();
		if (path.contains(v2_id))
		{
			if (firstV==v1_id || firstV==v2_id)
			{
				changed = true;
				if (firstV==v1_id || lastV==v2_id)

					// the whole path is inside the loop
					if ((firstV==v1_id && lastV==v1_id && path.size()==maxNumOfOccurrences*2+1) ||
							(firstV==v2_id && lastV==v2_id && path.size()==maxNumOfOccurrences*2-1) || 
							(firstV==v1_id && lastV==v2_id && path.size()==maxNumOfOccurrences*2) || 
							(firstV==v2_id && lastV==v1_id && path.size()==maxNumOfOccurrences*2) ) // all path is in the loop, but there is only one new path that matches
					{
						changed = updateSinglePathWithDoubleLoopNodes(path,v1_id,v2_id,newVerIDsV1,newVerIDsV2);
					}else
						path.clear();
				else
				{// only the start is inside the loop
					updatePathToRemoveLoopNodes(path, loopVs); 
					changed = true;
				}
			}else
			{ // start and ends outside the loop
				changed = updateSinglePathWithDoubleLoopNodes(path,v1_id,v2_id,newVerIDsV1,newVerIDsV2);
			}
		}
		if (changed)
			debugMes("\tpath changed from "+origPath+" to "+path,15);
		else
			debugMes("\tpath remains unchanged.", 15);
	}





	/**
	 * given this path, and the loop info, update the path to its single option.
	 * @param path
	 * @param v1_id
	 * @param v2_id
	 * @param newVerIDsV1
	 * @param newVerIDsV2
	 * @return
	 */
	private static boolean updateSinglePathWithDoubleLoopNodes(
			List<Integer> path, int v1_id, int v2_id, List<Integer> newVerIDsV1,
			List<Integer> newVerIDsV2) {
		boolean changed = false;
		for (int i=1 ; i<=path.size()-1 ; i++)
		{
			if (path.get(i).intValue()==v1_id) 
			{
				int j = newVerIDsV2.indexOf(path.get(i-1));
				if (j>=0)
				{
					path.set(i, newVerIDsV1.get(j+1));
					changed = true;
				}
			} else if (path.get(i).intValue()==v2_id) 
			{
				int j = newVerIDsV1.indexOf(path.get(i-1));
				if (j>=1)
				{
					path.set(i, newVerIDsV2.get(j));
					changed = true;
				}
			}
		}
		return changed;
	}

	/**
	 * print out the given error message, only if DEBUG=true
	 * @param mes Message
	 */
	private static void debugMes(String mes, int verbosityLevel)
	{
		//TODO: use general logging that can be leveraged across all classes.

		if (DEBUG && verbosityLevel<=BFLY_GLOBALS.VERBOSE_LEVEL)
		{
			if (USE_STDERR)
				System.err.println(mes);
			else if (LOG_STDERR)
				ERR_STREAM.println(mes);
		}

	}



	/**
	 * combine prefixes:
	 * calc for each v it's "depth" in terms of length of strings (from them on)
	 * draw all v's with the same depth
	 * sort on their set of parents
	 * draw all v's with same depth and same set of parents
	 * find subsets of those with same prefix
	 * create new node with prefix, connect accordingly.
	 * add the rest (those that removed the prefix) back into queue, with new depths
	 * @param graph
	 * @param compaction_round 
	 * @return
	 */
	
	/*
	private static void compactPrefixesBottomUp(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{

		debugMes("******   Compacting Prefixes BottomUp **********", 10);
		
		debugMes("setVerticesDepths()", 25);
		setVerticesDepths(graph);
		
		Comparator<Object> depthComparator = new SeqVertexDepthComparator();
		PriorityQueue<SeqVertex> dQueue = new PriorityQueue<SeqVertex>(graph.getVertexCount(),depthComparator );
		for (SeqVertex v : graph.getVertices())
			dQueue.add(v);

	
		SeqVertex v;
		ListComparator listComp = new ListComparator();
		
		boolean changed = false;
		
		for (int curD=0 ; curD<=MAX_DEPTH ; curD++) 
		{
			
			debugMes("COMPACTING_PREFIXES, curr_depth=" + curD, 20);
			
			// collect those sets of children that share exactly the same set of parents.
			TreeMap<List<SeqVertex>,Collection<SeqVertex>> curParents = new TreeMap<List<SeqVertex>,Collection<SeqVertex>>(listComp );
		

			while (!dQueue.isEmpty() && dQueue.peek().getDepth()==curD)
			{
				v = dQueue.poll();
				if (!graph.containsVertex(v))
					continue;

				List<SeqVertex> parents = getSortedParentList(graph,v); // list of parents of v, sorted by sequence
				
				if (!parents.isEmpty())
				{
					debugMes("curParents: "+curParents,20);
					if (!curParents.containsKey(parents))
					{
						debugMes(parents +" doesn't appear in curParents",20);
						curParents.put(parents,new HashSet<SeqVertex>());
					}
					debugMes("adding "+ v +" to "+curParents.get(parents),20);
					curParents.get(parents).add(v);
				}	

			}
			
			//look for groups of parents with identical children
			for (Collection<SeqVertex> parents : curParents.keySet())
			{
				// this collection has vertices with the same children
				Collection<SeqVertex> candidateNodes = curParents.get(parents);

				
				
				if (candidateNodes.size()==1)
					continue;

				
				debugMes("Collection of parents at same depth: " + candidateNodes, 20);
				
				// look for shared suffix
				boolean updateQueue = false;
				Collection<SeqVertex> updatedNodes = new HashSet<SeqVertex>();
				changed = compactPrefixRecursive(graph,candidateNodes,updatedNodes);

				for (SeqVertex ver : updatedNodes)
				{
					if (ver.getName().isEmpty())
					{
						debugMes("Need to update the queue. candidateNodes = "+updatedNodes,20);
						updateQueue = true;
					}
				}
				if (updateQueue)
				{
					for (SeqVertex ver : updatedNodes)
					{
						if (!ver.getName().isEmpty())
						{
							dQueue.add(ver);
							debugMes("adding "+ver+" to the queue, with depth "+ver.getDepth(),20);
						}
					}
				}
			}

		}

		
	}
	*/
	
	
	/*
	
	private static boolean mergeIdenticalPrefixesBottomUp(DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{

		debugMes("******   Merging Identical Prefixes BottomUp **********", 10);

		debugMes("setVerticesDepths()", 25);
		setVerticesDepths(graph);


		boolean changed = false;


		// group children by parents
		HashMap<List<SeqVertex>, List<SeqVertex>> shared_parents_to_children = new HashMap<List<SeqVertex>, List<SeqVertex>>(); 

		for (SeqVertex v : graph.getVertices()) {

			List<SeqVertex> parents = new ArrayList<SeqVertex>(graph.getPredecessors(v));

			if (! shared_parents_to_children.containsKey(parents)) {
				shared_parents_to_children.put(parents, new ArrayList<SeqVertex>());
			}
			shared_parents_to_children.get(parents).add(v);

		}

		for (List<SeqVertex> children_w_same_parents : shared_parents_to_children.values()) {

			if (children_w_same_parents.size() > 1) {

				Collection<SeqVertex> updatedNodes = new HashSet<SeqVertex>();

				changed = compactPrefixRecursive_v2(graph, children_w_same_parents, updatedNodes) || changed;

			}
		}

		return(changed);
	}
	
	*/



	
	/*
	private static boolean mergeIdenticalPrefixesBottomUp_v2 (DirectedSparseGraph<SeqVertex, SimpleEdge> graph)
	{

		// different from v1 version above in that all children of a given vertex are grouped for merge.
		// results is a very compact graph that still allows for generating proper sequences, and
		// disallowing for ties.  But, this graph is so compact that it generates larger number of 
		// repeat nodes, and slower unrolling times. 

		debugMes("******   Merging Identical Prefixes BottomUp **********", 10);

		debugMes("setVerticesDepths()", 25);
		setVerticesDepths(graph);


		boolean changed = false;

		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 20) {
			describeNodes(graph);
		}

		List<SeqVertex> all_vertices = new ArrayList<SeqVertex>(graph.getVertices());
		for (SeqVertex v : all_vertices) {

			if (! graph.containsVertex(v))
				continue;

			Collection<SeqVertex> all_children = graph.getSuccessors(v);

			// remove self vertex from list of children if we're in a direct loop:
			Collection<SeqVertex> children = new ArrayList<SeqVertex>();

			for (SeqVertex child : all_children) {
				if (! child.equals(v)) {
					children.add(child);
				}
			}


			if (children.size() > 1) {

				debugMes("## Targeting Merge Identical Prefix from parent: " + v 
						+ ", and children " + children, 20);

				Collection<SeqVertex> updatedNodes = new HashSet<SeqVertex>();
				changed = compactPrefixRecursive_v2(graph, children, updatedNodes) || changed;

			}
		}

		return(changed);
	}


	*/

	/**
	 * Given the graph, go over all vertices, and calculate their depth, as in distance from the roots 
	 * (maximal or minimal??) = doesn't matter as long as it's consistent. I chose maximal. 
	 * @param graph
	 */
	private static void setVerticesDepths(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		My_DFS dfs = new My_DFS(graph);
		dfs.runDFS2();
		
		List<SeqVertex> topBottom = getTopologicalOrder(graph, dfs);

		for (SeqVertex v : topBottom)
		{
			if (graph.inDegree(v)==0)
			{
				v.setDepth(0);
			}
			else
			{
				int d = -1;
				for (SeqVertex tv : graph.getPredecessors(v))
				{
					if (tv.getDepth() + tv.getNameKmerAdj().length() >d)
						d=tv.getDepth() + tv.getNameKmerAdj().length();
				}
				v.setDepth(d);
				
				
				// MAX_DEPTH global var being set here, used by prefix compaction method.
				if (d>MAX_DEPTH)
					MAX_DEPTH = d;
			}
		}

	}

	/**
	 * Given the graph, and the vertex v, return a sorted list of its parents
	 * @param graph
	 * @param v
	 * @return
	 */
	private static List<SeqVertex> getSortedParentList(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph, SeqVertex v) {
		
		
		List<SeqVertex> res = new ArrayList<SeqVertex>(graph.getPredecessors(v));
		SeqComparator verComp = new SeqComparator();
		Collections.sort(res, verComp);

		return res;
	}

	/**
	 * Given the graph, and the candidate nodes, look for shared prefixes of a single letter, 
	 * and move on.
	 * @param graph
	 * @param candidateNodes
	 * @param updateQueue 
	 */
	
	/*
	private static boolean compactPrefixRecursive(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Collection<SeqVertex> 	candidateNodes, 
			Collection<SeqVertex> updatedNodes) {
		
		
		debugMes("** compactPrefixRecursive: " + candidateNodes, 20);
		
		boolean changed = false;
		
		for (String l : LETTERS)
		{
			Collection<SeqVertex> vWithL = getVerticesWithFirstLetter(candidateNodes,l);
			if (vWithL.size()<=1)
				continue;

			debugMes("vWithL set based on l=" + l + ": " + candidateNodes, 20);
			

			// if there is a letter that has more than one vertex, create a new vertex with this letter
			changed = true;
			SeqVertex newV = new SeqVertex(getNextID(), l);
			
			// retain the identity of the vertices being collapsed here.
			newV.addIDsAsFirstPrevIDs(vWithL,LAST_REAL_ID);
			
			
			Collection<SeqVertex> new_vWithL = new HashSet<SeqVertex>(); 
			Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();

			for (SeqVertex v_in_vWithL : vWithL)
			{
				if (!graph.containsVertex(v_in_vWithL))
					continue;

				// create a new vertex with the first prevID as id
				SeqVertex newReplaceV_in_vWithL;
				if (!v_in_vWithL.getPrevVerIDs().isEmpty()    // node already subsumed another complex node
						&& 
						v_in_vWithL.getPrevVerIDs().firstElement().size()>1)
				{
					
					// just copying over the info into a new node 
					newReplaceV_in_vWithL = new SeqVertex(getNextID(), v_in_vWithL.getName());
					newReplaceV_in_vWithL.copyTheRest(v_in_vWithL);
					
					debugMes("compactPrefixRecursive/complex: Node: " + v_in_vWithL.getShortSeqWID() + " => " + newReplaceV_in_vWithL.getShortSeqWID(), 20);
					
				}
				else {
					
					newReplaceV_in_vWithL = v_in_vWithL.generateNewVerWithFirstIDasID(); 
					debugMes("compactPrefixRecursive/simple: Node: " + v_in_vWithL.getShortSeqWID() + " => " + newReplaceV_in_vWithL.getShortSeqWID(), 20);
				}
				
				// move all edges from and to the orig, to the new
				if (!newReplaceV_in_vWithL.equals(v_in_vWithL)) // they will be equal if the v_withL has no prevIDs, and only his original id
				{
					for (SimpleEdge e : graph.getOutEdges(v_in_vWithL))
					{
						removeEdges.add(e);
						graph.addEdge(new SimpleEdge(e.getWeight()), newReplaceV_in_vWithL, graph.getDest(e));
					}
					for (SimpleEdge e : graph.getInEdges(v_in_vWithL))
					{
						removeEdges.add(e);
						graph.addEdge(new SimpleEdge(e.getWeight()), graph.getSource(e), newReplaceV_in_vWithL);
					}
				}
				//replace it's location within vWithL
				new_vWithL.add(newReplaceV_in_vWithL);
			}

			for (SimpleEdge re : removeEdges)
			{
				debugMes("removing edge "+re+" between "+graph.getSource(re)+" and "+graph.getDest(re),20);
				graph.removeEdge(re);
			}

			for (SeqVertex rv : vWithL)
			{
				if (!new_vWithL.contains(rv))
				{
					debugMes("removing vertex (because new_vWithL doesn't contain it) "+rv,20);
					graph.removeVertex(rv);
				}
			}

			
			
			
			
			vWithL = new_vWithL;
			graph.addVertex(newV);
			debugMes("pulled the first letter from all vertices in "+vWithL+" to the new vertex "+newV,20);
			Vector<SeqVertex> removeVertices = new Vector<SeqVertex>();
			for (SeqVertex v1 : vWithL)
			{

				//				if (removeVertices.contains(v1) || !graph.containsVertex(v1))
				if (v1.isToBeDeleted() || !graph.containsVertex(v1))
					continue;
				removeEdges.clear();
				v1.increaseDepthByOne();

				////////////////////////////////////
				// reassign incoming edges to newV
				////////////////////////////////////
				
				
				for (SimpleEdge edgeToRemove : graph.getInEdges(v1))
				{
					double w2 = edgeToRemove.getWeight();
					SimpleEdge newE2 = null;
					SeqVertex v3 = graph.getSource(edgeToRemove);

					if (graph.findEdge(v3,newV)==null)
					{
						newE2 = new SimpleEdge(w2);
						graph.addEdge(newE2, v3,newV);  // edge reassignment
						debugMes("adding edge "+newE2+" between "+v3+" and "+newV,20);
					}else
					{
						newE2 = graph.findEdge(v3,newV);
						if (w2>newE2.getWeight())
						{
							//FIXME ?? do we want to add up the weights?
							debugMes("setting edge "+newE2+"'s weight from "+newE2.getWeight()+" to "+w2,20); 
							newE2.setWeight(w2);
						}
					}

					removeEdges.add(edgeToRemove);
					debugMes("removed edge "+edgeToRemove+" between "+graph.getSource(edgeToRemove)+" and "+graph.getDest(edgeToRemove),20);

				}
				
				///////////////////////////////////
				// handle outgoing edges (newE1)
				////////////////////////////////////
				
				
				
				if  (v1.getName().length()==1)
				{
					// single base successor node, just remove it and reassign out-edge.
					
					v1.removeFirstLetter();

					//go over all edges going out of v1, and move them to exit newV
					for (SeqVertex v0 : graph.getSuccessors(v1))
					{
						double w = graph.findEdge(v1,v0).getWeight();
						graph.addEdge(new SimpleEdge(w), newV,v0); // edge reassignments.
						debugMes("adding edge "+w+" between "+newV+" and "+v0,20);
					}
					debugMes("vertex "+v1+" is going to be removed",20);

					removeVertices.add(v1);
					v1.setToBeDeleted(true);

				}
				else if (v1.getName().length()<=KMER_SIZE && graph.outDegree(v1)==0) {
					
					// short terminal node less than a kmer size.
					// why do we need to handle this special use case?
					
					
					SeqVertex newV1; // needed only if this node is less than K in length
					
					v1.removeFirstLetter();
					
					Collection<SeqVertex> upV = graph.getPredecessors(v1);

					if (v1.getID()<=LAST_REAL_ID)
					{
						newV1 = new SeqVertex(getNextID(),v1.getName());
						graph.addVertex(newV1);
						removeVertices.add(v1);
						v1.setToBeDeleted(true);
					} else
						newV1 = v1; // what scenario?

					//go over all edges going into v1, and move them to exit newV
					if (upV.size()==1) // why only handling the case of a single parent?
					{
						for (SeqVertex upV1 : upV)  // only one here, no iteration needed
						{
							SimpleEdge oldE = graph.findEdge(upV1, v1);
							double w = oldE.getWeight();
							graph.addEdge(new SimpleEdge(w), newV,newV1);
							removeEdges.add(oldE); // already done above?
							debugMes("adding edge "+w+" between "+newV+" and "+newV1,20);
							debugMes("removing edge "+w+" between "+upV1+" and "+v1,20); 
							graph.addEdge(new SimpleEdge(1), v1, newV1);
						}
					}
				}
				else {
					
					// all other cases.
				
					double w = v1.removeFirstLetter();
					SimpleEdge newE1 = new SimpleEdge(w);
					graph.addEdge(newE1, newV,v1);
					debugMes("adding edge "+newE1+" between "+newV+" and "+v1,20);

				}

				for (SimpleEdge re : removeEdges)
				{
					graph.removeEdge(re);
				}
			}
			
			
			
			
			//try this out
			updatedNodes.clear();
			Set<SeqVertex> toAddTo_vWithL = new HashSet<SeqVertex>();
			
			
			// removing vertices targeted for deletion.
			for (SeqVertex rv : removeVertices)
			{
				graph.removeVertex(rv);
				debugMes("removed vertex "+rv,20);
				if (vWithL.contains(rv))
					vWithL.remove(rv);
				if (candidateNodes.contains(rv))
					candidateNodes.remove(rv);
			}

			// the restructuring to newV could result in new children available for further compaction.
			// check for other children of newV that are at the same depth and candidates for further compaction:
			for (SeqVertex vChild : graph.getSuccessors(newV))  
				if (!vWithL.contains(vChild) && vChild.hasAllSameParents(graph, vWithL))
						//vChild.getDepth()==curDepth)
					toAddTo_vWithL.add(vChild);
			
			
			for (SeqVertex vToAdd : toAddTo_vWithL)
				vWithL.add(vToAdd);

			for (SeqVertex vToAdd : vWithL)
			{
				updatedNodes.add(vToAdd);
			}


			if (vWithL.size()>1)
				compactPrefixRecursive(graph, vWithL, updatedNodes);

		}

		return changed;
	}

	
	*/
	
	
	
	/*
	private static boolean compactPrefixRecursive_v2(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			Collection<SeqVertex> 	candidateNodes, 
			Collection<SeqVertex> updatedNodes) {
		
		
		debugMes("** compactPrefixRecursive: " + candidateNodes, 20);
		
		boolean changed = false;
		
		for (String l : LETTERS)
		{
			Collection<SeqVertex> vWithL = getVerticesWithFirstLetter(candidateNodes,l);
			if (vWithL.size()<=1)
				continue;

			debugMes("vWithL set based on l=" + l + ": " + candidateNodes, 20);
			

			// if there is a letter that has more than one vertex, create a new vertex with this letter
			changed = true;
			SeqVertex newV = new SeqVertex(getNextID(), l);
			
			// retain the identity of the vertices being collapsed here.
			newV.addIDsAsFirstPrevIDs(vWithL,LAST_REAL_ID);
			
			
			
			// copy the current vertex list over to a new set of nodes.
			Collection<SeqVertex> new_vWithL = new HashSet<SeqVertex>(); 
			Vector<SimpleEdge> removeEdges = new Vector<SimpleEdge>();

			for (SeqVertex v_in_vWithL : vWithL)
			{
				if (!graph.containsVertex(v_in_vWithL))
					continue;

				// create a new vertex with the first prevID as id
				SeqVertex newReplaceV_in_vWithL;
				if (!v_in_vWithL.getPrevVerIDs().isEmpty()    // node already subsumed another complex node
						&& 
						v_in_vWithL.getPrevVerIDs().firstElement().size()>1)
				{
					
					// just copying over the info into a new node 
					newReplaceV_in_vWithL = new SeqVertex(getNextID(), v_in_vWithL.getName());
					newReplaceV_in_vWithL.copyTheRest(v_in_vWithL);
					
					debugMes("compactPrefixRecursive/complex: Node: " + v_in_vWithL.getShortSeqWID() + " => " + newReplaceV_in_vWithL.getShortSeqWID(), 20);
					
				}
				else {
					
					newReplaceV_in_vWithL = v_in_vWithL.generateNewVerWithFirstIDasID(); 
					debugMes("compactPrefixRecursive/simple: Node: " + v_in_vWithL.getShortSeqWID() + " => " + newReplaceV_in_vWithL.getShortSeqWID(), 20);
				}
				
				// move all edges from and to the orig, to the new
				if (!newReplaceV_in_vWithL.equals(v_in_vWithL)) // they will be equal if the v_withL has no prevIDs, and only his original id
				{
					for (SimpleEdge e : graph.getOutEdges(v_in_vWithL))
					{
						removeEdges.add(e);
						graph.addEdge(new SimpleEdge(e.getWeight()), newReplaceV_in_vWithL, graph.getDest(e));
					}
					for (SimpleEdge e : graph.getInEdges(v_in_vWithL))
					{
						removeEdges.add(e);
						graph.addEdge(new SimpleEdge(e.getWeight()), graph.getSource(e), newReplaceV_in_vWithL);
					}
				}
				//replace it's location within vWithL
				new_vWithL.add(newReplaceV_in_vWithL);
			}

			for (SimpleEdge re : removeEdges)
			{
				debugMes("removing edge "+re+" between "+graph.getSource(re)+" and "+graph.getDest(re),20);
				graph.removeEdge(re);
			}

			
			// remove the original vertices.
			for (SeqVertex rv : vWithL)
			{
				if (!new_vWithL.contains(rv))
				{
					debugMes("removing vertex (because new_vWithL doesn't contain it) "+rv,20);
					graph.removeVertex(rv);
				}
			}

			
			
			
			
			vWithL = new_vWithL;
			graph.addVertex(newV);
			debugMes("pulled the first letter from all vertices in "+vWithL+" to the new vertex "+newV,20);
			Vector<SeqVertex> removeVertices = new Vector<SeqVertex>();
			for (SeqVertex v1 : vWithL)
			{

				if (v1.isToBeDeleted() || !graph.containsVertex(v1))
					continue;
				removeEdges.clear();
				v1.increaseDepthByOne();

				////////////////////////////////////
				// reassign incoming edges to newV
				////////////////////////////////////
				
				
				for (SimpleEdge edgeToRemove : graph.getInEdges(v1))
				{
					double w2 = edgeToRemove.getWeight();
					SimpleEdge newE2 = null;
					SeqVertex v3 = graph.getSource(edgeToRemove);

					if (graph.findEdge(v3,newV)==null)
					{
						newE2 = new SimpleEdge(w2);
						graph.addEdge(newE2, v3,newV);  // edge reassignment
						debugMes("adding edge "+newE2+" between "+v3+" and "+newV,20);
					}else
					{
						newE2 = graph.findEdge(v3,newV);
						if (w2>newE2.getWeight())
						{
							//FIXME ?? do we want to add up the weights?
							debugMes("setting edge "+newE2+"'s weight from "+newE2.getWeight()+" to "+w2,20); 
							newE2.setWeight(w2);
						}
					}

					removeEdges.add(edgeToRemove);
					debugMes("removed edge "+edgeToRemove+" between "+graph.getSource(edgeToRemove)+" and "+graph.getDest(edgeToRemove),20);

				}
				
				///////////////////////////////////
				// handle outgoing edges (newE1)
				////////////////////////////////////
				
				
				
				if  (v1.getName().length()==1)
				{
					// single base successor node, just remove it and reassign out-edge.
					
					v1.removeFirstLetter();

					//go over all edges going out of v1, and move them to exit newV
					for (SeqVertex v0 : graph.getSuccessors(v1))
					{
						double w = graph.findEdge(v1,v0).getWeight();
						graph.addEdge(new SimpleEdge(w), newV,v0); // edge reassignments.
						debugMes("adding edge "+w+" between "+newV+" and "+v0,20);
					}
					debugMes("vertex "+v1+" is going to be removed",20);

					removeVertices.add(v1);
					v1.setToBeDeleted(true);

				}
				
				else {
					
					// all other cases.
				
					double w = v1.removeFirstLetter();
					SimpleEdge newE1 = new SimpleEdge(w);
					graph.addEdge(newE1, newV,v1);
					debugMes("adding edge "+newE1+" between "+newV+" and "+v1,20);

				}

				for (SimpleEdge re : removeEdges)
				{
					graph.removeEdge(re);
				}
			}
			
			
			
			
			//try this out
			updatedNodes.clear();
			Set<SeqVertex> toAddTo_vWithL = new HashSet<SeqVertex>();
			int curDepth = -1;
			// use this curDepth to decide if to add the children or not.
			if (!removeVertices.isEmpty())
				for (SeqVertex ver : vWithL)
					//					if (!removeVertices.contains(ver))
					if (!ver.isToBeDeleted())
						curDepth = ver.getDepth();
			
			
			// removing vertices targeted for deletion.
			for (SeqVertex rv : removeVertices)
			{
				graph.removeVertex(rv);
				debugMes("removed vertex "+rv,20);
				if (vWithL.contains(rv))
					vWithL.remove(rv);
				if (candidateNodes.contains(rv))
					candidateNodes.remove(rv);
			}

			// the restructuring to newV could result in new children available for further compaction.
			// check for other children of newV that are at the same depth and candidates for further compaction:
			for (SeqVertex vChild : graph.getSuccessors(newV))  {
				//debugMes("vChild: " + vChild+ ", vWithL: " + vWithL, 5);
				if (!vWithL.contains(vChild) && vChild.hasAllSameParents(graph, vWithL)) {
						//vChild.getDepth()==curDepth)
					toAddTo_vWithL.add(vChild);
				}
			}
			
			for (SeqVertex vToAdd : toAddTo_vWithL)
				vWithL.add(vToAdd);

			for (SeqVertex vToAdd : vWithL)
			{
				updatedNodes.add(vToAdd);
			}


			if (vWithL.size()>1)
				compactPrefixRecursive_v2(graph, vWithL, updatedNodes);

		}

		return changed;
	}

	*/
	
	
	
	/**
	 * Given the set of nodes, return a set of nodes that has the given letter l as a final letter
	 * @param candidateNodes
	 * @param l
	 * @return
	 */
	
	/*
	
	private static Collection<SeqVertex> getVerticesWithFirstLetter(
			Collection<SeqVertex> candidateNodes, String l) {
		Collection<SeqVertex> res = new HashSet<SeqVertex>();
		for (SeqVertex v : candidateNodes)
		{
			if (v.getName().startsWith(l))
				res.add(v);
		}
		return res;
	}

	*/

	// retrieve path list from first unshared node till the end (minus the final vertex)
	public static List<Integer> get_unshared_path_terminus(List<Integer> path_to_search, List<Integer> path_to_index) {

		debugMes("Path to search: " + path_to_search, 19);
		debugMes("Path to index: " + path_to_index, 19);

		Hashtable<Integer,Boolean> path_index = new Hashtable<Integer,Boolean>();
		for (Integer x : path_to_index) {
			path_index.put(x, new Boolean(true));
		}

		int unshared_path_pos = path_to_search.size(); // init to Infinity in essence, never reach this.
		for (int i = 0; i <= path_to_search.size()-2; i++) {
			if (! path_index.containsKey( path_to_search.get(i) ) ) {
				unshared_path_pos = i;
				break;
			}
		}

		List<Integer> unique_terminal_path = new Vector<Integer>();
		for (int i = unshared_path_pos; i <= path_to_search.size() -2; i++) {
			unique_terminal_path.add(path_to_search.get(i));
		}

		debugMes("Unique terminal path: " + unique_terminal_path, 19);

		return(unique_terminal_path);
	}	

	// see if any node is shared between the lists
	public static boolean paths_have_node_in_common (List<Integer> pathA, List<Integer> pathB) {


		Hashtable<Integer,Boolean> path_index = new Hashtable<Integer,Boolean>();
		for (int i = 0; i < pathA.size() - 1; i++) {
			path_index.put(pathA.get(i), new Boolean(true));
		}

		for (int i = 0; i < pathB.size() -1; i++) {
			if (path_index.containsKey( pathB.get(i))) {
				return(true);
			}
		}

		return(false);
	}	


	// see if any node other than the very last one is shared between the lists
	public static boolean paths_have_any_node_in_common (List<Integer> pathA, List<Integer> pathB, boolean include_sinks) {


		Hashtable<Integer,Boolean> path_index = new Hashtable<Integer,Boolean>();
		for (int i = 0; i < pathA.size() - 1; i++) {
			Integer node = pathA.get(i);
			if ( (! include_sinks) && node < 0) {
				continue; // sink node
			}
			path_index.put(node, new Boolean(true));
		}

		for (int i = 0; i < pathB.size() -1; i++) {
			Integer node = pathB.get(i);
			if (path_index.containsKey( node)) {
				// debugMes("Found node: " + node + " in common between paths: " + pathA + " and " + pathB, 10);
				return(true);
			}
		}

		return(false);
	}	



	public static String getPathMappingAsciiIllustration (
			final List<Integer> finalPath, 
			HashMap<PairPath,Integer> readPathsHashmap
			) {

		String ascii_illustration = "";

		for (int i = 0; i < finalPath.size(); i++) {
			ascii_illustration += "=";
		}
		ascii_illustration += "    PATH: " + finalPath + "\n";

		List<PairPath> readPaths = new ArrayList(readPathsHashmap.keySet());
		
		Collections.sort(readPaths, new Comparator<PairPath>() { // sort illustration by first node position in path
			public int compare(PairPath a, PairPath b) {
				Integer b_index = finalPath.indexOf(b.getFirstID());
				Integer a_index = finalPath.indexOf(a.getFirstID());
				return(a_index - b_index);
			}

		});

		for (PairPath read : readPaths) {

			char chars[] = new char[finalPath.size()];
			for (int i = 0; i < chars.length; i++) {
				chars[i] = ' ';
			}

			for (List<Integer> readPath : read.get_paths()) {
				for (Integer vertex_id : readPath) {
					int index = finalPath.indexOf(vertex_id);
					if (index >= 0) {
						chars[index] = '=';
					}
				}
			}
			for (int i = 0; i < chars.length; i++) {
				ascii_illustration += chars[i];
			}

			int read_counts = readPathsHashmap.get(read);
	

			ascii_illustration += "    Read: " + read.get_paths() + " read_support: " + read_counts + "\n";

		}


		return(ascii_illustration);

	}


	public static int count_pairpath_support(List<Integer> path, HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {

		HashMap<PairPath,Integer> pairPath_map = PathReads.get(path);

		int sum_reads = 0;

		for (PairPath p : pairPath_map.keySet()) {

			int read_count = pairPath_map.get(p);
			sum_reads += read_count;

		}

		return(sum_reads);	
	}


	public static HashMap<Integer, List<List<Integer>>> extractTripletsFromReads(HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash) {


		HashMap<Integer, List<List<Integer>>> tripletMapper = new HashMap<Integer, List<List<Integer>>>();

		for (Integer vertex_id : combinedReadHash.keySet()) {

			HashMap<PairPath,Integer> pp_hmap = combinedReadHash.get(vertex_id);

			for (PairPath pp : pp_hmap.keySet()) {

				
				
				List<List<Integer>> paths = pp.get_paths();
				for (List<Integer> read_path : paths) {
					if (read_path.size() < 3) {
						continue;
					}
					// iterate through triplets
					for (int i = 1; i < read_path.size()-1; i++) {
						Integer central_id = read_path.get(i);
						Integer left_id = read_path.get(i-1);
						Integer right_id = read_path.get(i+1);

						List<Integer> adjacency_path = new ArrayList<Integer>();
						adjacency_path.add(left_id);
						adjacency_path.add(central_id);
						adjacency_path.add(right_id);

						if (tripletMapper.containsKey(central_id)) {
							List<List<Integer>> triplet_list = tripletMapper.get(central_id);
							if (! triplet_list.contains(adjacency_path)) {
								triplet_list.add(adjacency_path);
								debugMes("Adding triplet adjacency_path to central node: " + central_id + " => " + adjacency_path, 17);
							}
							else {
								debugMes("triplet adjacency_path of node: " + central_id + " => " + adjacency_path + " already captured.", 17);
							}
						}
						else {
							List<List<Integer>> triplet_list = new ArrayList<List<Integer>>();
							triplet_list.add(adjacency_path);
							tripletMapper.put(central_id, triplet_list);
							debugMes("Setting initial triplet adjacency_path for central node: " + central_id + " => " + adjacency_path, 17);

						}



					}


				}

			}

		}
		
		
		return(tripletMapper);

	}

	
	public static HashMap<Integer, List<List<Integer>>> extractComplexPathPrefixesFromReads(HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash) {


		debugMes("-capturing path prefixes", 15);
		
		HashMap<Integer, List<List<Integer>>> nodeToComplexPathPrefix = new HashMap<Integer, List<List<Integer>>>();
		
		for (Integer vertex_id : combinedReadHash.keySet()) {

			HashMap<PairPath,Integer> pp_hmap = combinedReadHash.get(vertex_id);

			for (PairPath pp : pp_hmap.keySet()) {

				List<List<Integer>> paths = pp.get_paths();
				
				for (List<Integer> read_path : paths) {
					if (read_path.size() < 3) {
						continue;
					}
					// iterate through prefixes
					for (int i = read_path.size()-1; i >= 2; i--) {
						
						Integer node_id = read_path.get(i);	
						
						List<Integer> prefix_path = read_path.subList(0, i+1);
						
						
						if (! nodeToComplexPathPrefix.containsKey(node_id)) {
							nodeToComplexPathPrefix.put(node_id, new ArrayList<List<Integer>>());
						}
						
						if (! nodeToComplexPathPrefix.get(node_id).contains(prefix_path)) {
								nodeToComplexPathPrefix.get(node_id).add(prefix_path);
						}

					}


				}

			}

		}
		
		
		
		debugMes("-removing prefixes that are subpaths of other prefixes", 15);
		
		// remove paths that are subpaths
		for (Integer node_id : nodeToComplexPathPrefix.keySet()) {
			List<List<Integer>> prefixes = nodeToComplexPathPrefix.get(node_id);
			
			List<List<Integer>> prefixes_to_purge = new ArrayList<List<Integer>>();
			
			for (List<Integer> prefix : prefixes) {
				for (List<Integer> prefix2 : prefixes) {
					if (prefix != prefix2 && prefix2.size() > prefix.size() 
							&& Path.share_suffix_fully_contained(prefix, prefix2)) {
						prefixes_to_purge.add(prefix);
					}
				}
			}
			
			for (List<Integer> prefix : prefixes_to_purge) {
				prefixes.remove(prefix);
			}
		}
		

		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 17) {
			for (List<List<Integer>> extendedTripletPathsList : nodeToComplexPathPrefix.values()) {

				for (List<Integer> extendedTripletPath : extendedTripletPathsList) {
					debugMes("EXTENDED_TRIPLET_CAPTURED: " + extendedTripletPath, 17);
				}
			}
		}
		
		
		return(nodeToComplexPathPrefix);

	}
	
	
	
	public static Boolean tripletSupported(List<List<Integer>> triplet_list, List<Integer> triplet) {

		for (List<Integer> t_list : triplet_list) {
			debugMes("Checking triplet list: " + t_list + " comparing to query triplet: " + triplet, 15);
			if (t_list.get(0).equals(triplet.get(0))
					&& 
					t_list.get(1).equals(triplet.get(1))
					&&
					t_list.get(2).equals(triplet.get(2))
					) {

				return(true);
			}	
		}
		return(false);


	}


	public static List<Integer> ensure_path_has_sinks(List<Integer> path) {

		List<Integer> new_path = new ArrayList<Integer>(path);

		if (new_path.get(0) != -1) {
			new_path.add(0, -1);
		}
		if (new_path.get(new_path.size()-1) != -2) {
			new_path.add(-2);
		}

		return(new_path);


	}


	public static HashMap<List<Integer>, Pair<Integer>> reduce_cdhit_like (
			HashMap<List<Integer>, Pair<Integer>> FinalPaths_all, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<List<Integer>,HashMap<PairPath,Integer>> PathReads) {



		////////////////////////////////////////////////////////////////////////////
		debugMes("\n\n****  CD-HIT style path collapsing at end of run.\n\n", 15);


		Vector<FinalPaths> path_vec = new Vector<FinalPaths>();

		DecimalFormat df = new DecimalFormat("#.#");
		for (List<Integer> path : FinalPaths_all.keySet())
		{
			String seq = getPathSeq(graph,path);

			FinalPaths f = new FinalPaths(path, seq);
			path_vec.add(f);
		}
		

		MAX_SEQ_LEN_DP_ALIGN  = ALL_VS_ALL_MAX_DP_LEN; // temporarily replace

		Collections.sort(path_vec); // sort paths by length of sequence descendingly

		// examine sequence CD-HIT -style, remove those that lack sufficient variation
		HashMap<FinalPaths,Boolean> filtered = new HashMap<FinalPaths,Boolean>();


		List<List<Integer>> removeSimilarPaths =  new ArrayList<List<Integer>>(); 

		for (int i = 0; i < path_vec.size()-1; i++) {

			if (filtered.containsKey(path_vec.get(i))) {
				// path filtered, cannot use it as evidence for filtering smaller sequences.
				continue;
			}

			List<Integer> path_i = path_vec.get(i).path;

			List<Integer> path_i_w_sinks = ensure_path_has_sinks(path_i);

			for (int j = i + 1; j < path_vec.size(); j++) {


				if (filtered.containsKey(path_vec.get(j))) {
					// path filtered, cannot use it as evidence for filtering smaller sequences.
					continue;
				}

				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
					System.err.print("\r[" + i + "," + j + "]   ");
				}
				
				List<Integer> path_j = path_vec.get(j).path;
				List<Integer> path_j_w_sinks = ensure_path_has_sinks(path_j);

				
				
				
				String seq_i = path_vec.get(i).sequence;
				String seq_j = path_vec.get(j).sequence;

				/*
				int index1 = seq_i.length()-1;
				int index2 = seq_j.length()-1;
				debugMes("ALL-VS-ALL: (" + i + "," + j + " of " + path_vec.size() + ") checking for similarity the two paths: "+path_i+ 
						"(len="+seq_i.length()+");"+path_j+"(len="+seq_j.length()+")",10);

				 */


				//if (finalSeqsAreTooSimilar(seq_i, seq_j)) {

				if (twoPathsAreTooSimilar(graph, path_i_w_sinks, path_j_w_sinks)) {
					
					debugMes("\n\n*** REDUCE: they are TOO SIMILAR! ***\n\n",10);	

					int rIndex = removeTheLesserSupportedPath(seq_i, seq_j, path_i, path_j, removeSimilarPaths, PathReads);

					//int rIndex = removeTheShorterPath(path1S,path2S,path1,path2,removeSimilarPaths,PathReads);

					if (rIndex == 1) {// the first path was removed

						filtered.put(path_vec.get(i), true);
						debugMes("\tRemoving (" + i + ") seq in pair", 18);
					}
					else {
						filtered.put(path_vec.get(j), true);
						debugMes("\tRemoving (" + j + ") second seq in pair", 18);
					}

				}
				else 
					debugMes("\n\n*** REDUCE: they are PLENTY DIFFERENT ***\n\n", 15);
			}
		}


		for (FinalPaths path2Remove : filtered.keySet())
		{
			debugMes("REDUCE-STAGE: The final path "+path2Remove+" was removed because it was too close to another path",10);
			FinalPaths_all.remove(path2Remove.path);

		}

		ALL_VS_ALL_MAX_DP_LEN = MAX_SEQ_LEN_DP_ALIGN; // back to original setting

		return(FinalPaths_all);

	}

	public static HashMap<List<Integer>, Pair<Integer>> remove_identical_subseqs (
			HashMap<List<Integer>, Pair<Integer>> FinalPaths_all, 
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			HashMap<List<Integer>,HashMap<PairPath,Integer>> PathReads) {



		
		Vector<FinalPaths> path_vec = new Vector<FinalPaths>();

		DecimalFormat df = new DecimalFormat("#.#");
		int count = 0;
		for (List<Integer> path : FinalPaths_all.keySet())
		{
			count++;
			debugMes("-reconstructing sequence for path[: " + count + " of " + FinalPaths_all.keySet().size() + "]: " + path, 15);
			String seq = getPathSeq(graph,path);

			FinalPaths f = new FinalPaths(path, seq);
			path_vec.add(f);
		}


		////////////////////////////////////////////////////////////////////////////
		debugMes("\n\n****  Removing identical subsequences among: " + path_vec.size() + " paths.\n\n", 10);

		
		Collections.sort(path_vec); // sort paths by length of sequence descendingly

		// examine sequence CD-HIT -style, remove those that lack sufficient variation
		HashMap<FinalPaths,Boolean> filtered = new HashMap<FinalPaths,Boolean>();


		List<List<Integer>> removeSimilarPaths =  new ArrayList<List<Integer>>(); 

		for (int i = 0; i < path_vec.size()-1; i++) {

			if (filtered.containsKey(path_vec.get(i))) {
				// path filtered, cannot use it as evidence for filtering smaller sequences.
				continue;
			}

			List<Integer> path_i = path_vec.get(i).path;

			List<Integer> path_i_w_sinks = ensure_path_has_sinks(path_i);

			for (int j = i + 1; j < path_vec.size(); j++) {


				if (filtered.containsKey(path_vec.get(j))) {
					// path filtered, cannot use it as evidence for filtering smaller sequences.
					continue;
				}

				List<Integer> path_j = path_vec.get(j).path;
				List<Integer> path_j_w_sinks = ensure_path_has_sinks(path_j);

				String seq_i = path_vec.get(i).sequence;
				String seq_j = path_vec.get(j).sequence;

				int index1 = seq_i.length()-1;
				int index2 = seq_j.length()-1;
				
				if (BFLY_GLOBALS.VERBOSE_LEVEL >= 15) {
					System.err.print("\r[" + i + "," + j + "]    ");
				}
				else {
					debugMes("ALL-VS-ALL: (" + i + "," + j + " of " + path_vec.size() + ") checking for identical subseqs between the two paths: "+path_i+ 
							"(len="+seq_i.length()+");"+path_j+"(len="+seq_j.length()+")",16);
				}

				if (seq_i.indexOf(seq_j) >= 0) { 
					filtered.put(path_vec.get(j), true);
					debugMes("\t** Removing (" + j + ") seq in pair, contains " + i, 15);
				}
				else if (seq_j.indexOf(seq_i) >= 0) {
					filtered.put(path_vec.get(i), true);
					debugMes("\t** Removing (" + i + ") seq in pair, contains " + j, 15);
				}

			}

		}


		for (FinalPaths path2Remove : filtered.keySet())
		{
			debugMes("REDUCE-STAGE: The final path "+path2Remove+" was removed because it was too close to another path",10);
			FinalPaths_all.remove(path2Remove.path);

		}


		return(FinalPaths_all);

	}


	public static Boolean finalSeqsAreTooSimilar (String seq_i, String seq_j) {

		// note, to check for perfectly identical sequence clusters, could run cd-hit like so:
		// cd-hit-est -o cdhit -c 1 -i comp.allProbPaths.fasta -p 1 -d 0 -b 1


		if ( (seq_i.length() > MAX_SEQ_LEN_DP_ALIGN && seq_j.length() > MAX_SEQ_LEN_DP_ALIGN) 
				||
				seq_i.length() > 100000   // problems can arise in the alignment code if either seq is longer
				||
				seq_j.length() > 100000) {


			// zipper:  Just get rid of those that are truly nearly identical.

			AlignmentStats stats = ZipperAlignment.doZipperAlignment("A", seq_i, "B", seq_j);
			int mismatches = stats.mismatches;
			debugMes("-zipper reports: " + mismatches + " mismatches between seqs.", 18);
			if (mismatches <= 2)
				return(true);
			else
				return(false);
		}	
		else {

			Alignment alignment;

			if (SMITH_WATERMAN_ALIGN_FLAG) {

				debugMes("-running Smith-Waterman alignment of path sequences", 15);
				alignment = NWalign.run_SW_alignment("A", seq_i, "B", seq_j, 4, -5, 10, 1);


			}
			else {
				// Needleman Wunsch Global Alignment is default
				debugMes("-running Needleman-Wunsch alignment of path sequences", 15);
				alignment = NWalign.run_NW_alignment("A", seq_i, "B", seq_j, 4, -5, 10, 1);   //NW locks up or takes too long with very long sequences (eg. 40kb align to 6kb)


			}


			int max_diffs_in_window = AlignmentStats.get_max_diffs_in_window(alignment, DIFFS_WINDOW_SIZE);
			debugMes (new jaligner.formats.Pair().format(alignment), 10);
			debugMes("Max diffs found in alignment window: " + max_diffs_in_window, 10);


			if (max_diffs_in_window <= MAX_FINAL_DIFFS_IN_WINDOW) {
				return (true);
			}


		}

		return(false);

	}			

	private static String get_pathName_string (List<Integer> path,
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {


		/*
		int startI = 0, endI = path.size(); 
		if (path.get(0)== ROOT.getID())
			startI++;
		if (path.get(path.size()-1)== T_VERTEX.getID())
			endI--;

		String pathName;
		if (MISO_OUTPUT) {
			pathName = "[";
			int iSeqL=0,j=0;
			for (int vi=startI; vi<endI; vi++){
				iSeqL = getSeqVertex(graph, path.get(vi)).getName().length();
				pathName = pathName + path.get(vi)+":"+j+"-"+(j+iSeqL-1);
				if (vi<endI-1)
					pathName = pathName.concat(" ");
				j+=iSeqL;
			}
			pathName = pathName +"]";
		} else
			pathName = ""+path.subList(startI, endI);

		 */

		int startI = 0, endI = path.size(); 
		if (path.get(0)== ROOT.getID())
			startI++;
		if (path.get(path.size()-1)== T_VERTEX.getID())
			endI--;


		String pathName,degenString="";
		SeqVertex v;
		
		
		if (MISO_OUTPUT) {

			DecimalFormat f0 = new DecimalFormat("#");

			pathName = "[";
			int iSeqL=0,j=0;
			for (int vi=startI; vi<endI; vi++){
				v = getSeqVertex(graph, path.get(vi));
				iSeqL = v.getName().length();
				if (vi != startI){
					iSeqL -= (KMER_SIZE-1);
				}
				
				//int node_id = v.getOrigButterflyID();
				int node_id = v.getID();
				String path_node = "" + node_id;
				
				/*
				if (xStructuresResolvedByTriplets.containsKey(node_id) && ! xStructuresResolvedByTriplets.get(node_id)) {
					path_node = "@" + path_node + "@!";
				}
				*/
				
				pathName = pathName + path_node +":"+j+"-"+(j+iSeqL-1);

				if (vi<endI-1)
					pathName = pathName.concat(" ");
				if (v.getDegenerativeFreq().size()>0)
				{
					for (int di=0; di<v.getDegenerativeFreq().size(); di++)
					{
						degenString = degenString + "{("+(j+v.getDegenerativeLocations().elementAt(di))+")";
						degenString = degenString + v.getDegenerativeLetters().elementAt(di).elementAt(0)+":"+
								f0.format(v.getDegenerativeFreq().elementAt(di).elementAt(0))+" ";
						degenString = degenString + v.getDegenerativeLetters().elementAt(di).elementAt(1)+":"+
								f0.format(v.getDegenerativeFreq().elementAt(di).elementAt(1)) +"}";
					}
				}
				j+=iSeqL;

			}
			pathName = pathName +"]";
		} else
			pathName = ""+path.subList(startI, endI);


		if (USE_DEGENERATE_CODE) {
			pathName = pathName + " SNPs="+degenString;
		}

		pathName += " " + path;
		
		
		return(pathName);
	}

	
	 public static HashMap<Integer,HashMap<PairPath,Integer>> getComponentReads(Integer component_id, Set<SeqVertex> comp, 
			 	HashMap<Integer,HashMap<PairPath,Integer>> combinedReadHash) {
		 
		 HashSet<Integer> node_ids_in_component = new HashSet<Integer>();
		 
		 for (SeqVertex v : comp) {
			 node_ids_in_component.add(v.getID());
		 }
		
		 List<Integer> component_node_list = new ArrayList<Integer>(node_ids_in_component);
		 
		 HashMap<Integer,HashMap<PairPath,Integer>> componentReadHash = new HashMap<Integer,HashMap<PairPath,Integer>>();
		 
		 for (Integer node_id : combinedReadHash.keySet()) {
			 HashMap<PairPath,Integer> pp_map = combinedReadHash.get(node_id);
			 for (PairPath p : pp_map.keySet()) {
				 if (p.haveAnyNodeInCommon(component_node_list)) {
					 if (! componentReadHash.containsKey(node_id)) {
						 componentReadHash.put(node_id, new HashMap<PairPath,Integer>());
					 }
					 
					 componentReadHash.get(node_id).put(p, pp_map.get(p));
					 debugMes("Subcomponent: " + component_id + ", adding pairpath: " + p, 15);
				 }
			 }
		 }
		 
		 return(componentReadHash);
	 }
	
	
	 /* bad idea... Cannot rely on node depth positions as they're imperfectly ordered
	 
	 public static List<Integer> constructSinglePathFromPairPathList (List<PairPath> pairpath_list,
			 final DirectedSparseGraph<SeqVertex, SimpleEdge> graph, 
				DijkstraDistance<SeqVertex,SimpleEdge> dijkstraDis
			 ) {
		 
		 
		 
		 
		 HashSet<Integer> extracted = new HashSet<Integer>();
		 ArrayList<SeqVertex> extractedVertices = new ArrayList<SeqVertex>();

		 for (PairPath pp : pairpath_list) {

			 extracted.addAll(pp.getPath1());

			 if(pp.hasSecondPath())
				 extracted.addAll(pp.getPath2());
		 }

		 for(Integer id : extracted)
		 {
			 extractedVertices.add(getSeqVertex(graph, id));
		 }
		 //extractedVerticesIDs.addAll(extracted);

		 // sort according to topological order of BTFL graph 
		 Collections.sort(extractedVertices, new SeqVertexFinishTimeComparator());

		 String node_id_list_text = "";
		 for (SeqVertex v : extractedVertices) {
			 node_id_list_text += v.getID() + " ";
		 }
		 debugMes("Extracted sorted vertices: " + node_id_list_text + "\n", 10);

		 DijkstraShortestPath dsp = new DijkstraShortestPath(graph);

		 // Fill in any gaps
		 // path=[L_1]; For each i in 1:length(node_list)
		 int j = 0;
		 int num_vertices = extractedVertices.size(); // note, extractedVertices grows in size during iterations below.
		 while(j < num_vertices - 1)
		 {		

			 //System.out.println(j);
			 SeqVertex current = extractedVertices.get(j);
			 SeqVertex next = extractedVertices.get(j + 1);

			 // 3.3.1 if L_i == L_(i+1) then nothing
			 // -There are no duplicates since extractedVertices was created from
			 //  building the set of extracted vertex IDs

			 // 3.3.2 else if exists an edge from L_i to L_(i+1) then nothing(?)
			 if(graph.getSuccessors(current).contains(next)) {
				 j++;
				 continue;
			 }

			 // 3.3.3 else find a single path (p = L_i,..., L_(i+1)): 
			 // append all P_j (j=2:end) to our path

			 //List<SimpleEdge> sp = org.jgrapht.alg.DijkstraShortestPath.findPathBetween((Graph)graph, current, next);

			 List<SimpleEdge> sp = dsp.getPath(current, next);
			 debugMes("Found shortest path between " + current.getID() + " and " + next.getID() + ":", 10);
			 ArrayList<SeqVertex> toAdd = new ArrayList<SeqVertex>();
			 for(SimpleEdge edge : sp) {
				 SeqVertex v = graph.getDest(edge);
				 toAdd.add(v);
				 debugMes("\t" + v.getID(), 10);
			 }
			 toAdd.remove(next);
			 extractedVertices.addAll(toAdd);

			 j++;
		 }

		 ArrayList<Integer> extractedVerticesIDs = new ArrayList<Integer>();
		 Collections.sort(extractedVertices, new SeqVertexFinishTimeComparator());
		 for(SeqVertex v: extractedVertices)
		 {
			 //System.out.println("Adding vertex with ID: " + v.getID());
			 extractedVerticesIDs.add(v.getID());
		 }


		 
		 return(extractedVerticesIDs);
		 
	 }
	 */
	 
	 
	 public static boolean twoPairPathsAreTooFarAwayInGraph (PairPath pp_i, PairPath pp_j, DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		 
		 SeqVertex last_vertex_i = getSeqVertex(graph, pp_i.getLastID());
		 SeqVertex first_vertex_j = getSeqVertex(graph, pp_j.getFirstID());
		 
		 
		 int last_vertex_i_start_time = last_vertex_i._node_depth;
		 int first_vertex_j_start_time = first_vertex_j._node_depth;
		
		 
		 /*
		 System.err.println("Pairpaths: " + pp_i + " vs. " + pp_j);
		 System.err.println("Start times: nodes [" + pp_i.getLastID() + ","  + pp_j.getFirstID() + "], times: " + i_start_time + " vs. " + j_start_time);
		
		 */
		 
		 if (first_vertex_j_start_time - last_vertex_i_start_time > MAX_VERTEX_DISCOVERY_TIME_DIFF_ALLOW_COMPARE) 
			 return(true);
		 else 
			 return(false);
		 
		 
	 }
	 
	 
	 
	 
	 
	 
}

// End TransAssembly.java
//
//		



