
import java.util.*;

import edu.uci.ics.jung.graph.DirectedSparseGraph;


public class Path {

	private List<Integer> vertex_id_listing;
	
	private int pathNodeID;
	private static int pathNodeCounter = 0;
	
	private static HashMap<String,Path> _pathNodeIndexer = new HashMap<String,Path>();
	
	int _tmp_score;
	
	
	public Path (List<Integer> path) {
		this.vertex_id_listing = path;
		pathNodeCounter++;
		pathNodeID = pathNodeCounter;
		
		_pathNodeIndexer.put(getPathNodeID(), this);
	}
	
	public List<Integer> get_vertex_list () {
		return(this.vertex_id_listing);
	}
	
	public int size() {
		return (this.vertex_id_listing.size());
	}
	
	public String getPathNodeID() {
		return("PN" + pathNodeID);
	}
	
	public String toString () {
		return(this.getPathNodeID() + "::" + this.vertex_id_listing);
	}
	
	
	
	///////////////////////////
	// various static methods
	///////////////////////////
	
	public static Path retrievePathNodeObj (String pathNodeID) {
		return(_pathNodeIndexer.get(pathNodeID));
	}
	
	
	public static Boolean contains_node_id (List<Integer> path, Integer node_id) {
		
		if (path.contains(node_id)) {
			return(true);
		}
		else {
			return(false);
		}
	}
	
	
	public static Boolean contains_any_node_id(List<Integer> path, Collection<Integer> node_ids) {
		
		for (Integer node_id : node_ids)  {
			
			if (contains_node_id(path, node_id)) {
				return(true);
				
			}
			
		}
		
		return(false); //none found.
	}


	public static int countNumNodesNotUnique(List<Integer> reconstructed_path) {
		
		int non_unique_count = 0;
		
		HashMap<Integer,Integer> counter = new HashMap<Integer,Integer>();
		for (Integer node_id : reconstructed_path)
		{
			if (! counter.containsKey(node_id)) {
				counter.put(node_id, 1);
			}
			else {
				counter.put(node_id, counter.get(node_id)+1);
			}
		}
		
		for (Integer node_id : reconstructed_path)
		{
			if (counter.get(node_id) > 1)
				non_unique_count++;
		}
		
		return(non_unique_count);
	}


	public static int countNumOutOfOrder(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			List<Integer> path) {
	
		int count_out_of_order = 0;
		
		for (int i = 1; i < path.size(); i++) {
			SeqVertex prev_v = SeqVertex.retrieveSeqVertexByID(path.get(i-1));
			SeqVertex curr_v = SeqVertex.retrieveSeqVertexByID(path.get(i));
			
			if (prev_v.getNodeDepth() > curr_v.getNodeDepth()) {
				// out of order
				count_out_of_order++;
			}
		}
		
		return(count_out_of_order);
		
	}


	public static int countNumOrigNodesNotUnique(List<Integer> path) {
		
		List<Integer> orig_id_path = new ArrayList<Integer>();
		for (Integer node_id : path) {
			SeqVertex v = SeqVertex.retrieveSeqVertexByID(node_id);
			orig_id_path.add(v.getOrigButterflyID());
		}
		
		int non_unique_count = 0;
		
		HashMap<Integer,Integer> counter = new HashMap<Integer,Integer>();
		for (Integer node_id : orig_id_path)
		{
			if (! counter.containsKey(node_id)) {
				counter.put(node_id, 1);
			}
			else {
				counter.put(node_id, counter.get(node_id)+1);
			}
		}
		
		for (Integer node_id : orig_id_path)
		{
			if (counter.get(node_id) > 1)
				non_unique_count++;
		}
		
		return(non_unique_count);
		
	}
	
	
	public static HashMap<Integer, Integer> getRepeatNodesAndCounts(List<Integer> path) {

		// want max node occurrence per path (don't count it twice if it shows up in two pair paths... harder to solve that one.
		HashMap<Integer,Integer> max_node_ids_n_counts = new HashMap<Integer,Integer>();


		HashMap<Integer,Integer> node_ids_n_counts = new HashMap<Integer,Integer>();
		for (Integer node_id : path) {
			if (node_ids_n_counts.containsKey(node_id)) {
				node_ids_n_counts.put(node_id,  node_ids_n_counts.get(node_id)+1);
			}
			else {
				node_ids_n_counts.put(node_id, 1);
			}
		}

		// set to max
		for (Integer node_id : node_ids_n_counts.keySet()) {
			Integer count = node_ids_n_counts.get(node_id);
			if (count > 1) {
				max_node_ids_n_counts.put(node_id, count);
			}
		}
		
		return(max_node_ids_n_counts);
		
	}


	public static double getUnrolledEdgeWeightSum(
			DirectedSparseGraph<SeqVertex, SimpleEdge> graph,
			List<Integer> reconstructed_path) {
		
		double unrolled_weight_sum = 0;
		
		HashSet<SimpleEdge> seen_edges = new HashSet<SimpleEdge>();
		
		for (int i = 1; i < reconstructed_path.size(); i++ ) {
			
			SeqVertex before_v = SeqVertex.retrieveSeqVertexByID(reconstructed_path.get(i-1));
			SeqVertex after_v = SeqVertex.retrieveSeqVertexByID(reconstructed_path.get(i));
			
			SimpleEdge se = graph.findEdge(before_v, after_v);
			
			if (se == null)
				throw new RuntimeException("Error, no edge found between " + before_v + " and " + after_v);
			
			if (! seen_edges.contains(se)) {
			
				// don't count the edge more than once!!!  
				
				unrolled_weight_sum += se.get_repeat_unroll_weight();
				
				seen_edges.add(se);
			}
		}
		
		
		return(unrolled_weight_sum);
		
	}


	public static boolean hasTerminalNonSelfRepeat(List<Integer> path) {
		
		if (path.size() <= 1) {
			return(false);
		}
		
		if (path.indexOf(path.get(path.size()-1)) != path.size()-1) {
			return(true);
		}
		else {
			return(false);
		}
	}

	
	public static List<Integer> create_path_from_single_node_id (Integer node_id) {
		
		List<Integer> path = new ArrayList<Integer>();
		path.add(node_id);
		
		return(path);
		
		
	}
	
	public static List<List<Integer>> create_list_of_paths_from_single_node_id(Integer node_id) {
		List<Integer> path = create_path_from_single_node_id(node_id);
		
		List<List<Integer>> path_list = new ArrayList<List<Integer>>();
		path_list.add(path);
		
		return(path_list);
	}
	
	public static List<List<Integer>> create_list_of_paths_including_path(List<Integer> path) {
		
		List<List<Integer>> path_list = new ArrayList<List<Integer>>();
		path_list.add(path);
		
		return(path_list);
	}
	
	public static List<List<Integer>> prepend_node_id_to_paths (Integer node_id, List<List<Integer>> path_list) {
		
		for (List<Integer> path : path_list) {
			path.add(0, node_id);
		}
		
		return(path_list);
		
	}
	
	
	public static List<List<Integer>> clone (List<List<Integer>> path_list) {
		
		List<List<Integer>> cloned_list = new ArrayList<List<Integer>>();
		
		for (List<Integer> path : path_list) {
			
			List<Integer> cloned_path = new ArrayList<Integer>(path);
			cloned_list.add(cloned_path);
			
		}
		
		return(cloned_list);
		
	}
	
	
	
	
	public static List<List<Integer>> create_empty_path_list () {
		List<List<Integer>> path_list = new ArrayList<List<Integer>>();
		
		path_list.add(new ArrayList<Integer>());
		
		return(path_list);
	}


	public static List<Integer> collapse_compatible_paths(
			List<List<Integer>> chain_i_path_list) {
		
		//System.err.println("Collapsing paths: "  + chain_i_path_list);
		
		List<List<Integer>>  uncollapsed_paths = Path.clone(chain_i_path_list);
		
		List<Integer> complete_path = uncollapsed_paths.remove(0);
		
		while (! uncollapsed_paths.isEmpty()) {
			
			List<List<Integer>> tmp_uncollapsed_paths = new ArrayList<List<Integer>>();
			
			boolean merged_flag = false;
			
			for (List<Integer> path : uncollapsed_paths) {
				if (PairPath.haveAnyNodeInCommon(complete_path, path)) {
					if (! PairPath.individual_paths_are_compatible(complete_path, path)) {
						throw new RuntimeException("Error, paths share node but are not compatible: " + complete_path + " and " + path);
					}
					
					complete_path = Path.mergeCompatiblePaths(complete_path, path);
					
					merged_flag = true;
				}
				else {
					tmp_uncollapsed_paths.add(path);
				}
			}
			
			if (! merged_flag) {
				throw new RuntimeException("Error, no merged paths in this round, so_far_collapsed_path: " + complete_path 
						+ ", remaining unmerged paths: " + tmp_uncollapsed_paths);
			}
			
			
			uncollapsed_paths.clear();
			if (! tmp_uncollapsed_paths.isEmpty()) {
				uncollapsed_paths.addAll(tmp_uncollapsed_paths);
			}
			
		}
		
		return(complete_path);
		
		
	}


	private static List<Integer> mergeCompatiblePaths(
			List<Integer> pathA, List<Integer> pathB) {
		
		//System.err.println("merging paths: " + pathA + " and " + pathB);
		
		if (pathA.equals(pathB)) {
			return(pathA);
		}
		else if (Path.pathA_contains_pathB(pathA, pathB)) {
			return(pathA);
		}
		else if (Path.pathA_contains_pathB(pathB, pathA)) {
			return(pathB);
		}
		else if (pathA.indexOf(pathB.get(0)) > 0) {
			// extend complete path
			return(Path.extend_pathA_by_pathB(pathA, pathB));
				
		}
		else if (pathB.indexOf(pathA.get(0)) > 0) {
			return(Path.extend_pathA_by_pathB(pathB, pathA));
		
		}
		else {
			throw new RuntimeException("can't merge supposedly compatible paths: " 
					+ pathA + " and " + pathB);
		}
	}


	private static List<Integer> extend_pathA_by_pathB(List<Integer> pathA,
			List<Integer> pathB) {
		
		
		int pathB_start_index = pathA.indexOf(pathB.get(0));
		if (pathB_start_index < 0) {
			throw new RuntimeException("Error, cannot extend path, as pathB doesn't start within A: A=" + pathA + ", B=" + pathB);
		}
		
		List<Integer> extended_path = new ArrayList<Integer>();
		for (int i = 0; i < pathB_start_index; i++) {
			extended_path.add(pathA.get(i));
		}
		
		extended_path.addAll(pathB);
		
		return(extended_path);
	}


	private static boolean pathA_contains_pathB(List<Integer> pathA,
			List<Integer> pathB) {
		
		int pathB_start_index = pathA.indexOf(pathB.get(0));
		if (pathB_start_index < 0) {
			return(false);
		}
		
		int i,j;
		
		for (i = pathB_start_index,  j = 0; i < pathA.size() && j < pathB.size(); i++, j++) {
			
			if ( ! pathB.get(j).equals(pathA.get(i)) ) {
				return(false);
			}
		}
		
		if (j != pathB.size()) {
			return(false); // whole path not included.
		}
		
		// passed all tests. A contains B 
		return(true);
		
	}


	public static List<Integer> collapse_compatible_pair_paths(List<PairPath> paths) {
		
		// note this is over-simplified, as the pp are separated into their individual paths
		// but methods using this should have been smarter about deciding what to collapse based on pp constraints.
		
		List<List<Integer>> all_individual_paths = new ArrayList<List<Integer>>();
		
		for (PairPath pp : paths) {
			all_individual_paths.add(pp.getPath1());
			if (pp.hasSecondPath()) {
				all_individual_paths.add(pp.getPath2());
			}
		}
		
		return(Path.collapse_compatible_paths(all_individual_paths));
	}


	public static List<List<Integer>> collapse_compatible_paths_to_min_set(
			List<List<Integer>> chain_i_path_list) {


		//System.err.println("Collapsing paths: "  + chain_i_path_list);

		List<List<Integer>>  uncollapsed_paths = Path.clone(chain_i_path_list);

		List<List<Integer>> complete_paths = new ArrayList<List<Integer>>();
		complete_paths.add(uncollapsed_paths.remove(0));

		while (! uncollapsed_paths.isEmpty()) {

			List<Integer> complete_path = complete_paths.get(complete_paths.size()-1);

			List<List<Integer>> tmp_uncollapsed_paths = new ArrayList<List<Integer>>();

			boolean merged_flag = false;

			for (List<Integer> path : uncollapsed_paths) {
				if (PairPath.haveAnyNodeInCommon(complete_path, path)
						&&
						PairPath.individual_paths_are_compatible(complete_path, path)
						) {


					complete_path = Path.mergeCompatiblePaths(complete_path, path);

					merged_flag = true;
				}
				else {
					tmp_uncollapsed_paths.add(path);
				}
			}

			if (! merged_flag) {
				// seed a new entry
				complete_paths.add(tmp_uncollapsed_paths.remove(0));
			}


			uncollapsed_paths.clear();
			if (! tmp_uncollapsed_paths.isEmpty()) {
				uncollapsed_paths.addAll(tmp_uncollapsed_paths);
			}

		}

		return(complete_paths);
		
	}


	public static boolean pathA_contains_pathB_allowRepeats(List<Integer> pathA, List<Integer> pathB) {


		for (int i = 0; i < pathA.size(); i++) {
			if (pathB.get(0).equals(pathA.get(i))) {
				
				int matches = 0;
				
				// have a first position match.  See if it extends to a complete match:
				for (int pathB_pos=0, pathA_pos = i; 
						pathA_pos < pathA.size() 
						 && pathB_pos < pathB.size();
						pathA_pos++, pathB_pos++) {
				
					Integer pathA_id = pathA.get(pathA_pos);
					Integer pathB_id = pathB.get(pathB_pos);
					
					
					if (pathB_id.equals(pathA_id)) {
						matches++;
						
						/*
						System.err.println("" + pathB_pos + "=" + pathB_id
								 + " vs. " + pathA_pos + "=" + pathA_id + " matches.");
								 
						 */
					}
					else {
						
						/*
						System.err.println("" + pathB_pos + "=" + pathB_id
								 + " vs. " + pathA_pos + "=" + pathA_id + " MISmatches.");
						*/
						
						break;
					}
				}
				
				//System.err.println("Got match count: " + matches + " vs. path size " + pathB.size());
				
				if (matches == pathB.size()) {
					// got it.
					//System.err.println("** found match OK");
					return(true);
					
				}
				
			}
		}
		
		//System.err.println("*!* no match.");		
		return(false); // not found.
		
	}

	public static PathOverlap pathB_extends_pathA_allowRepeats(List<Integer> pathB, List<Integer> pathA, HashSet<Integer> repeat_node_ids) {

		int best_match_count = -1;
		
		PathOverlap path_overlap = new PathOverlap(0,0);
		
		for (int i = pathA.size()-1; i >= pathA.size() - pathB.size() && i >= 0; i--) {
			if (pathB.get(0).equals(pathA.get(i))) {

				int matches = 0;
				int overlap_len = 0;


				// have a first position match.  See if it extends to a complete match:
				for (int pathB_pos=0, pathA_pos = i; 
						pathA_pos < pathA.size() 
						&& pathB_pos < pathB.size();
						pathA_pos++, pathB_pos++) {

					Integer pathA_id = pathA.get(pathA_pos);
					Integer pathB_id = pathB.get(pathB_pos);

					overlap_len++;
					if (pathB_id.equals(pathA_id)) {
						if (! repeat_node_ids.contains(pathA_id)) {
							// dont count it if it's a repeat node.
							matches++;
						}
						/*
						System.err.println("" + pathB_pos + "=" + pathB_id
								+ " vs. " + pathA_pos + "=" + pathA_id + " matches.");
						*/
					}
					else {

						/*
						System.err.println("" + pathB_pos + "=" + pathB_id
								+ " vs. " + pathA_pos + "=" + pathA_id + " MISmatches.");
						*/

						matches = -1;
						break;
					}
				}

				//System.err.println("Got match count: " + matches + " vs. path size " + pathB.size());

				if (matches > best_match_count) {
					best_match_count = matches;
					path_overlap.match_length = overlap_len;
					path_overlap.match_score = matches;
					
				}

			}
		}

		return(path_overlap);

	}

	public static boolean share_suffix_fully_contained(List<Integer> path_i,
			List<Integer> path_j) {
		
		
		for (int i = path_i.size() - 1, j = path_j.size() - 1; 
				i >= 0 && j >= 0; 
				i--, j--) {
			
			if (path_i.get(i) != path_j.get(j)) {
				return(false);
			}
			
		}
		
		return(true);		
		
	}

}
