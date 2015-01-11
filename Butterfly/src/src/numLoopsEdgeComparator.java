import java.util.Comparator;

import edu.uci.ics.jung.graph.DirectedSparseGraph;


public class numLoopsEdgeComparator implements Comparator<Object> {

	
	DirectedSparseGraph<SeqVertex, SimpleEdge> _graph;
	
	public numLoopsEdgeComparator (DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		
		_graph = graph;
		
	}
	
	
	@Override
	public int compare(Object o1, Object o2) {
		SimpleEdge e1 = (SimpleEdge)o1;
		SimpleEdge e2 = (SimpleEdge)o2;
		
		int l1 = e1.getNumOfLoopsInvolved();
		int l2 = e2.getNumOfLoopsInvolved();
		
		double s1 = e1.getWeight();
		double s2 = e2.getWeight();
		
		
		// prioritize by number of loops involved in.
		if( l1 < l2 )
			return 1;
		else if( l1 > l2 )
			return -1;
		else
		{
			// prioritize by lower edge weight
			
			if (s1 > s2)
				return 1;
			else if (s1 < s2)
				return -1;
			else {
				
				
				return(0);
				
				/*
			
				//FIXME: use hierarchical layout to set node depths.
			
				// prioritize by absolute difference in node depths (order descending in abs(delta_depth)
				SeqVertex v1_a = _graph.getSource(e1);
				SeqVertex v1_b = _graph.getDest(e1);
				int abs_delta_e1_depth = Math.abs(v1_a.getDepth() - v1_b.getDepth());
				
				SeqVertex v2_a = _graph.getSource(e2);
				SeqVertex v2_b = _graph.getDest(e2);
				int abs_delta_e2_depth = Math.abs(v2_a.getDepth() - v2_b.getDepth());
				
				if (abs_delta_e1_depth < abs_delta_e2_depth) {
					return(1);
				}
				else if (abs_delta_e1_depth > abs_delta_e2_depth) {
					return(-1);
				}
				else {
					return(0);
				}
				*/
				
				
			}
		}
		
	}

}
