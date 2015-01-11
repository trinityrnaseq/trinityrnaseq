import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

public class PathReadSupportComparator implements Comparator<Object> {

	HashMap<List<Integer>, HashMap<PairPath, Integer>> _path_reads;
	
	HashMap<List<Integer>, Integer> _read_support;
	
	public PathReadSupportComparator (HashMap<List<Integer>, HashMap<PairPath, Integer>> PathReads) {
		_path_reads = PathReads;
		_read_support = new HashMap<List<Integer>, Integer>();
	}
	
	
	@Override
	public int compare(Object o1, Object o2) {
		
		List<Integer> path1 = (List<Integer>) o1;
		List<Integer> path2 = (List<Integer>) o2;
		
		int read_support_path1 = get_read_support(path1);
		int read_support_path2 = get_read_support(path2);
		
		if (read_support_path1 < read_support_path2) {
			return(-1);
		}
		else if (read_support_path1 > read_support_path2) {
			return(1);
		}
		else {
			return(0);
		}
	
	}
	
	
	private int get_read_support (List<Integer> path) {
		
		if (_read_support.containsKey(path)) {
			int read_count = _read_support.get(path);
			return(read_count);
		}
		else {
			// sum up all pair paths from reads 
			HashMap<PairPath,Integer> pairPath_map = _path_reads.get(path);
			
			int sum_reads = 0;
			
			for (PairPath p : pairPath_map.keySet()) {
				
				int read_count = pairPath_map.get(p);
				sum_reads += read_count;
				
			}
		
			_read_support.put(path, sum_reads);
			return(sum_reads);
		
		}
		
		
	}

	
	
	
}
