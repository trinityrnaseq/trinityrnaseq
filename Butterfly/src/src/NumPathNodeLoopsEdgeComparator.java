import java.util.Comparator;


public class NumPathNodeLoopsEdgeComparator implements Comparator<SimplePathNodeEdge> {

	@Override
	public int compare(SimplePathNodeEdge e1, SimplePathNodeEdge e2) {
		
		if (e1._num_loops < e2._num_loops) {
			return(1);
		}
		else if (e1._num_loops > e2._num_loops) {
			return(-1);
		}
		else {
			
			// same number of loops, prioritize by smallest weight
			if (e1.weight < e2.weight) {
				return(-1);
			}
			else if (e1.weight > e2.weight) {
				return(1);
			}
			else {
				return(0);
			}
		}
		
	}

}
