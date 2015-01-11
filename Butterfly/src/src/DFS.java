import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import edu.uci.ics.jung.graph.DirectedSparseGraph;


public class DFS
{
	private Map<SeqVertex,String> _colors;
	public final static String WHITE = "white";
	public final static String BLACK = "black";
	public final static String GRAY = "gray";
	protected DirectedSparseGraph<SeqVertex,SimpleEdge> _graph;
	protected Set<SeqVertex> _roots;
	/**
	 * Constructor for the DFS object
	 * @param roots 
	 */
	public DFS(DirectedSparseGraph<SeqVertex,SimpleEdge> graph, Set<SeqVertex> roots) 
	{ 
		_colors = new HashMap<SeqVertex,String>();
		_graph = graph;
		_roots = roots;
		for (SeqVertex v : roots)
		{
			_colors.put(v,WHITE);
		}
	}

	/**
	 * Gets the color attribute of the DFS object
	 */
	public String getColor(SeqVertex v)
	{
		String color = _colors.get(v);
		return color != null ? color : WHITE;
	}

	/**
	 * Visit an edge
	 */
	private void visitEdge(SimpleEdge e,
			Visitor visitor)
	{
		visitor.discoverEdge(e);

		SeqVertex v = _graph.getDest(e);

		if (getColor(v) == WHITE)
		{
			visitVertex(v, visitor);
		}

		visitor.finishEdge(e);
	}

	/**
	 * Visit a vertex:
	 * color it gray, and explore all its descendants
	 */
	private void visitVertex(SeqVertex v,
			Visitor visitor)
	{
		_colors.put(v, GRAY);

		visitor.discoverVertex(v); // mark the discovery time of this vertex

		Iterator<SimpleEdge> edges = _graph.getOutEdges(v).iterator();
		while (edges.hasNext())
		{
			SimpleEdge e = (SimpleEdge) edges.next();
			visitEdge(e, visitor);
		}

		visitor.finishVertex(v); // mark the finishing time of this vertex

		_colors.put(v, BLACK);
	}

//	/**
//	 * visit - Visits the graph
//	 */
//	public void visit(SeqVertex root,
//			Visitor visitor)
//	{
//		_colors.clear();
//		visitor.discoverGraph(_graph);
//
//		visitVertex(root, visitor);
//
//		visitor.finishGraph(_graph);
//	}

	/**
	 * visit - Visits all nodes in the graph.
	 * go over all roots, if they are white, visit this vertex further (in recursion)
	 */
	public void visit(Visitor visitor ) {
		//    _colors.clear();
		visitor.discoverGraph( _graph );

		Iterator<SeqVertex> vertices = _roots.iterator();
		while (vertices.hasNext()) {
			SeqVertex v = (SeqVertex) vertices.next();

			if (_colors.get( v ) == WHITE) {
				visitVertex(v, visitor );
			}
		}

		visitor.finishGraph( _graph );
	}


}



