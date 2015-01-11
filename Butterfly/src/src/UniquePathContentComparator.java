import java.util.*;


public class UniquePathContentComparator implements Comparator<List<Integer>> {

	
	HashMap<List<Integer>, Integer> path_unique_read_content;
	HashMap<List<Integer>, Integer> seqLengthMap;
	
	public UniquePathContentComparator(
			List<List<Integer>> paths,
			HashMap<PairPath, Boolean> pp_used,
			HashMap<List<Integer>, HashMap<PairPath, Integer>> finalPathsToContainedReads, HashMap<List<Integer>, Integer> seqLengthMap) {
	
		this.seqLengthMap = seqLengthMap;
		
		path_unique_read_content = new HashMap<List<Integer>,Integer>();
		
		for (List<Integer> path : paths) {
			
			path_unique_read_content.put(path, 0);
			
			if (finalPathsToContainedReads.containsKey(path)) {
				HashMap<PairPath,Integer> compatible_pp_hmap = finalPathsToContainedReads.get(path);

				for (PairPath pp : compatible_pp_hmap.keySet()) {
					if (! pp_used.containsKey(pp)) {
						int pp_read_support = compatible_pp_hmap.get(pp);
						path_unique_read_content.put(path, path_unique_read_content.get(path) + pp_read_support);
					}
				}
				
			}
			
		}
		
		
		
	}

	public int compare (List<Integer> path1, List<Integer> path2) {
		
		int unique_read_content_path1 = this.path_unique_read_content.get(path1);
		int unique_read_content_path2 = this.path_unique_read_content.get(path2);
		
		if (unique_read_content_path1 < unique_read_content_path2) {
			return(-1);
		}
		else if (unique_read_content_path1 > unique_read_content_path2) {
			return(1);
		}
		else {
			// tied for read support.
			// sort by sequence length
			
			int seq_length_path1 = this.seqLengthMap.get(path1);
			int seq_length_path2 = this.seqLengthMap.get(path2);

			if (seq_length_path1 < seq_length_path2) {
				return(-1);
			}
			else if (seq_length_path1 > seq_length_path2) {
				return(1);
			}
			else {
				return(0);
			}
		}
	}

	public List<List<Integer>> remove_paths_without_unique_read_content(
			List<List<Integer>> paths) {
		
		List<List<Integer>> paths_retained = new ArrayList<List<Integer>>();
		
		for (List<Integer> path : paths) {
			if (path_unique_read_content.get(path) > 0) {
				paths_retained.add(path);
			}
		}
		
		return(paths_retained);
		
	}

	public int unique_count(List<Integer> path) {
		
		return(this.path_unique_read_content.get(path));
		
	}
	
	
	
	
	
}
