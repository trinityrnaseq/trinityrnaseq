import edu.uci.ics.jung.graph.DirectedSparseGraph;

import java.util.*;

public class TopologicalSort {

	

	
	/*
	 * 
	 * http://en.wikipedia.org/wiki/Topological_sorting
	 * Kahn, 1962 method:
	 * 
	 * 
	L ← Empty list that will contain the sorted elements
	S ← Set of all nodes with no incoming edges
	while S is non-empty do
	    remove a node n from S
	    add n to tail of L
	    for each node m with an edge e from n to m do
	        remove edge e from the graph
	        if m has no other incoming edges then
	            insert m into S
	if graph has edges then
	    return error (graph has at least one cycle)
	else 
	    return L (a topologically sorted order)
	
	*/
	
	
	
	public static List<SeqVertex> topoSortSeqVerticesDAG(
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph) {
		
		//L ← Empty list that will contain the sorted elements
		List<SeqVertex> L = new ArrayList<SeqVertex>();
		
		// S ← Set of all nodes with no incoming edges
		List<SeqVertex> S = new ArrayList<SeqVertex>();
		
		for (SeqVertex v : seqvertex_graph.getVertices()) {
			if (seqvertex_graph.getPredecessorCount(v) == 0) {
				S.add(v);
			}
			
			
		}
		
		
		/*
		 * while S is non-empty do
	    remove a node n from S
	    add n to tail of L
	    for each node m with an edge e from n to m do
	        remove edge e from the graph
	        if m has no other incoming edges then
	            insert m into S
		 */
		
		HashSet<SimpleEdge> edges_to_ignore = new HashSet<SimpleEdge>();
		
		
		
		// while S is non-empty do
		while (! S.isEmpty()) {
			// remove a node n from S
			SeqVertex n = S.remove(0);
			
			// add n to tail of L
			L.add(n);
			
			//for each node m with an edge e from n to m do
			for (SeqVertex m : seqvertex_graph.getSuccessors(n)) {
				SimpleEdge se = seqvertex_graph.findEdge(n, m);
				
				// remove edge e from the graph
				edges_to_ignore.add(se);
			
			
				// (determine num of parents still actively connected to m)
				int num_remaining_connected_preds = get_num_remaining_connected_predecessors(seqvertex_graph, m, edges_to_ignore);
				// if m has no other incoming edges then
				if (num_remaining_connected_preds == 0) {
					//insert m into S
					S.add(m);
				}
			}

			if (S.isEmpty()) {
				for (SeqVertex v : seqvertex_graph.getVertices()) {
					if ( (! L.contains(v)) ) {

						int pred_count = seqvertex_graph.getPredecessorCount(v);
						System.err.println("seq vertex " + v + " not selected yet and has pred count: " + pred_count);
					}

				}

			}
			
			
		}
		
		if (! edges_to_ignore.containsAll(seqvertex_graph.getEdges())) {
			for (SimpleEdge e : seqvertex_graph.getEdges()) {
				if (! edges_to_ignore.contains(e)) {
					System.err.println("ERROR: after topo sort, still have edge unaccounted for: " + e);
				}
			}
			throw new RuntimeException("Error, graph contains at least one cycle and is not a DAG!");
		}

		// assign node depths
		int depth = -1;
		for (SeqVertex v : L) {
			depth++;
			v.setDepth(depth);
			v.setNodeDepth(depth);
		}

		
		return(L);
		
	}

	private static int get_num_remaining_connected_predecessors(
			DirectedSparseGraph<SeqVertex, SimpleEdge> seqvertex_graph,
			SeqVertex m, HashSet<SimpleEdge> edges_to_ignore) {
		
		
		int num_remaining_connected_preds = 0;
		for (SeqVertex p : seqvertex_graph.getPredecessors(m)) {
			SimpleEdge se2 = seqvertex_graph.findEdge(p, m);
			if (! edges_to_ignore.contains(se2)) {
				num_remaining_connected_preds++;
			}
		}
		
		
		return(num_remaining_connected_preds);
		
	}
	
	
	
	
	
}
