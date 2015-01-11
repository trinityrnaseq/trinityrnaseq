import java.util.*;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.DirectedSparseGraph;


public class My_DFS
{
	private Map<SeqVertex,Integer> _colors;
	public final static int WHITE = 0;
	public final static int BLACK = 1;
	public final static int GRAY = 2;
	public final static int GREEN = 3;
	public final static int DARKGRAY = 4;
	
	protected DirectedSparseGraph<SeqVertex,SimpleEdge> _graph;
	protected Set<SeqVertex> _roots;
	
	protected int _time;
	protected Map<SeqVertex,Number> _discovery;
	protected Map<SeqVertex,Number> _finishing;
	
	DijkstraShortestPath<SeqVertex, SimpleEdge> _dp;
	
	boolean local_debug = false;

	/**
	 * Constructor for the DFS object
	 * Initialize all vertices to WHITE
	 * @param roots 
	 */
	public My_DFS(DirectedSparseGraph<SeqVertex,SimpleEdge> graph) 
	{ 
		_graph = graph;
		initDFS();
	}

	/**
	 * init all variables.
	 */
	public void initDFS()
	{
		_colors = new HashMap<SeqVertex,Integer>();
		for (SeqVertex v : _graph.getVertices()) {
			_colors.put(v,WHITE);
			v._node_depth = -1;
			if (local_debug) {
				System.err.println("[dfs init]: " + v.getShortSeqWconnectingIDs(_graph));
			}
		}
		
		_time=0;
		_discovery = new HashMap<SeqVertex,Number>();
		_finishing = new HashMap<SeqVertex,Number>();
		_dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(_graph);
	}
	
	
	/*
	
	public void runDFS()
	{
		initDFS();
		for (SeqVertex v : _graph.getVertices()) {
			if (_graph.inDegree(v)==0 && getColor(v)==WHITE) {
				visitVertex(v);
			}
		}
			
	}
	
	*/
	
	
	public void runDFS2()
	{
		initDFS();
		for (SeqVertex v : _graph.getVertices()) {
			if (_graph.inDegree(v)==0 && getColor(v)==WHITE) {
				if (local_debug) {
					System.err.println("INITIATING runDFS2() at: " + v.getShortSeqWconnectingIDs(_graph));
				}
				visitVertex2(v, 0);
			}
		}
		
		
		//TODO: find the minimum depth
		Integer min_depth = null;
		
		for (SeqVertex v : _graph.getVertices()) {
			if (min_depth == null)
				min_depth = v.getNodeDepth();
			else 
				min_depth = Math.min(min_depth, v.getNodeDepth());
		}
		
		//TODO: adjust all depth values so that min depth starts at zero
		
		for (SeqVertex v : _graph.getVertices()) {
			int node_depth = v.getNodeDepth();
			int adj_node_depth = node_depth - min_depth;
			v.setDepth(adj_node_depth);
			v.setNodeDepth(adj_node_depth);
			
		}
		
		//refine_depths_for_branched_parents_and_linear_below();
		
		//update_start_node_depths_bottom_up();
		
		employ_regular_down_only_DFS_search();
		
		
	}
	
	
	
	private void employ_regular_down_only_DFS_search() {
		
		//System.err.println("employ_regular_down_only_DFS_search()");
		
		List<SeqVertex> vertices = new ArrayList(_graph.getVertices());
		
		// sort by current depth settings.
		Collections.sort(vertices, new Comparator<SeqVertex>() {
			
			public int compare(SeqVertex v1, SeqVertex v2) {
				if (v1.getNodeDepth() < v2.getNodeDepth()) {
					return(-1);
				}
				else if (v1.getNodeDepth() > v2.getNodeDepth()) {
					return(1);
				}
				else {
					return(0);
				}
			}
			
		});
		
		
		// reinit colors:
		_colors = new HashMap<SeqVertex,Integer>();
		for (SeqVertex v : vertices) {
			v.setNodeDepth(-1);
			v.setDepth(-1);
			setColor(v,WHITE);
		
		}
		


		for (SeqVertex v : vertices) {
			if (getColor(v) == WHITE) {
				int depth = 0;
				if (_graph.getPredecessorCount(v) != 0) {
					depth = get_depth_based_on_parents(v);
				}	
				v.setDepth(depth);
				v.setNodeDepth(depth);
				setColor(v,GRAY);
				for (SeqVertex s : _graph.getSuccessors(v)) {
					down_only_visit_vertex(s, 1);
				}
				setColor(v, BLACK);

			}


			
		}
		
	}

	private int get_depth_based_on_parents(SeqVertex v) {
		int depth = -1;
		for (SeqVertex p : _graph.getPredecessors(v)) {
			depth = Math.max(depth, p.getNodeDepth());
		}
		
		depth++; // one more than parents.
		
		return(depth);
	}
	
	private void down_only_visit_vertex(SeqVertex v, int depth_to_assign) {
		if (getColor(v) != WHITE) {
			// already visited
			return;
		}
		if (has_white_parent(v)) {
			// delay
			return;
		}
		
		setColor(v, GRAY);
		v.setDepth(depth_to_assign);
		v.setNodeDepth(depth_to_assign);  // should unify these two depth vars.
		
		for (SeqVertex s : _graph.getSuccessors(v)) {
			down_only_visit_vertex(s, depth_to_assign+1);
		}
		
		setColor(v, BLACK); // done visiting vertex.
		
	}

	private boolean has_white_parent(SeqVertex v) {
		
		for (SeqVertex p : _graph.getPredecessors(v)) {
			if (getColor(p) == WHITE) {
				return(true);
			}
		}
		
		return(false);
	}

	private void refine_depths_for_branched_parents_and_linear_below() {
		
		List<SeqVertex> vertices = new ArrayList(_graph.getVertices());
		
		// sort by current depth settings.
		Collections.sort(vertices, new Comparator<SeqVertex>() {
			
			public int compare(SeqVertex v1, SeqVertex v2) {
				if (v1.getNodeDepth() < v2.getNodeDepth()) {
					return(-1);
				}
				else if (v1.getNodeDepth() > v2.getNodeDepth()) {
					return(1);
				}
				else {
					return(0);
				}
			}
			
		});
		
		
		
		for (SeqVertex v : vertices) {
			
			if (_graph.getPredecessorCount(v) > 1) {
				
				int max_parent_depth = v.get_max_parent_node_depth(_graph);
				
				if (v.getNodeDepth() != max_parent_depth + 1) {
					int node_depth = max_parent_depth + 1;
					
					//System.err.println("-resetting up-branched node " + v + " depth to " + node_depth);
					
					v.setNodeDepth(node_depth);
					
					
					
					// trickle down
					SeqVertex child_v = v;
					while (_graph.getSuccessorCount(child_v) == 1) {
						child_v = _graph.getSuccessors(child_v).iterator().next();
						if (_graph.getPredecessorCount(child_v) > 1) {
							break; // branched upwards, not linear stretch
						}
						node_depth++;
						
						//System.err.println("-trickle down: resetting " + child_v +  " depth to: " + node_depth);
						
						child_v.setNodeDepth(node_depth);
						
						
					}
					
				}
			}
			
		}
		
		
	}

	public void update_start_node_depths_bottom_up() {
		
		for (SeqVertex v : _graph.getVertices()) {
			if (_graph.inDegree(v) == 0 && _graph.outDegree(v) == 1) {
				// start node that is itself unbranched.
				// walk it until we find a branch point.
				
				List<SeqVertex> node_list = new ArrayList<SeqVertex>();
				node_list.add(v);
				
				SeqVertex next_v = v;
				
				while (_graph.outDegree(next_v) == 1 && _graph.inDegree(next_v) == 1) {
					node_list.add(next_v);
					next_v = _graph.getSuccessors(next_v).iterator().next();
				}
				
				if (_graph.getSuccessors(next_v).iterator().hasNext()) { // next is a branched (before or after) node.
					int branched_node_depth = _graph.getSuccessors(next_v).iterator().next()._node_depth;
					Collections.reverse(node_list);
					for (SeqVertex u : node_list) {
						branched_node_depth--;
						u._node_depth = branched_node_depth;
					}
				}
				
			}
		}
		
	}
	
	
	
	
	
//	private void updateCircularity() {
//		DijkstraShortestPath<SeqVertex, SimpleEdge> dp = new DijkstraShortestPath<SeqVertex, SimpleEdge>(_graph);
//
//		for (SeqVertex v : _graph.getVertices())
//			if (v.isInCircle())
//			{
//				assert(dp.getDistance(v, v)!=null); //v is reachable from itself
//				for (SimpleEdge e : dp.getPath(v, v))
//				{
//					e.setInCircle(true);
//					_graph.getDest(e).setInCircle(true);
//				}
//			}
//	}

	
	
	
	/**
	 * Visit a vertex:
	 * color it gray, and explore all its descendants
	 */
	
	/*
	private void visitVertex(SeqVertex v)
	{
		setColor(v, GRAY);
		_time++;
		// mark the discovery time of this vertex
		_discovery.put(v, _time);
		v.setDFS_DiscoveryTime(_time);
		
		//System.err.println("DFS: visiting node: " + v.getID());
		
		// determine node depth based on max (parent node depth) + 1
		Integer node_depth = null;
		Collection<SeqVertex> predecessors = _graph.getPredecessors(v);
		for (SeqVertex p : predecessors) {
			
			
			if (BFLY_GLOBALS.VERBOSE_LEVEL >= 25) {
				System.err.println("(dfs) Node " + v.getShortSeqWID() +
						" has predecessor " + p.getShortSeqWID());
			}
			
			if (node_depth == null)
				node_depth = p.getNodeDepth();
			else 
				node_depth = Math.max(node_depth, p.getNodeDepth());
				
		}
		if (node_depth == null) // has no predecessors
			node_depth = 0;
		else
			node_depth += 1;
		
		v._node_depth = node_depth;
		if (BFLY_GLOBALS.VERBOSE_LEVEL >= 25) {
			System.err.println("DFS: setting node: " + v.getID() + " depth: " + node_depth + " => " + v.getShortSeqWID());
		}
		
		for (SeqVertex u : _graph.getSuccessors(v))
		{
			if (getColor(u)==WHITE)
				visitVertex(u);
			else if (getColor(u)==GRAY){ // we have reached a circle.
				u.setInCircle(true);
				assert(_dp.getDistance(u, v)!=null); //v is reachable from itself
				if (v.equals(u)) //self loop
					_graph.findEdge(v, v).setInCircle(true);
				else //other loops
					for (SimpleEdge e : _dp.getPath(u, v))
					{
						e.setInCircle(true);
						_graph.getDest(e).setInCircle(true);
					}

			}
		}
		
		setColor(v, BLACK);
		_time++;
		// mark the finishing time of this vertex
		_finishing.put(v, _time);
		v.setDFS_FinishingTime(_time);
		
	
	}
	*/
	
	
	
	/**
	 * Visit a vertex:
	 * color it gray, and explore all its descendants, then the unvisited ancestors
	 */
	private void visitVertex2(SeqVertex v, int recommended_depth)
	{
		
		
		if (getColor(v) != WHITE)  { // already visited.
			if (local_debug) {
				System.err.println("\talready visited...: " + v.getShortSeqWconnectingIDs(_graph) + " having color: " + getColor(v));
			}
			return;
		}
		
		
		//////////////////////////////////////////////
		//  Setting depth based on recommendation
		//////////////////////////////////////////////
		
		v.setDepth(recommended_depth);
		v.setNodeDepth(recommended_depth);
		
		if (local_debug) {
			System.err.println("DFS ENTRY for: " + v.getShortSeqWconnectingIDs(_graph));
		}

		setColor(v, GRAY); // initial visit, depth set based on recommendation
		_time++;
		// mark the discovery time of this vertex
		_discovery.put(v, _time);
		v.setDFS_DiscoveryTime(_time);
		if (local_debug) {
			System.err.println("\tset discovery time for V" + v.getID()  + " to " + "T" + _time);
		}
		
		///////////////////
		// Visit children
		///////////////////
		
		
		for (SeqVertex s : _graph.getSuccessors(v)) {
			if (local_debug) {
				System.err.println("\tvisiting successor: " + s.getID() + " of " + v.getID());
			}
			visitVertex2(s, recommended_depth + 1);
		}
		
		//////////////////////////
		//  Visit parent nodes:
		/////////////////////////
		
		for (SeqVertex p : _graph.getPredecessors(v)) {
			
			if (local_debug) {
				System.err.println("\tvisiting predecessor: " + p.getID() + " of " + v.getID());
			}
			visitVertex2(p, recommended_depth - 1);
		}

		///////////////////////////////////////////////////////
		// back from having visited all parents
		// redefine depth based on max of parent depth value.
		///////////////////////////////////////////////////////
		
		
		/*
		Integer max_up_depth = null;
		
		Collection<SeqVertex> predecessors = _graph.getPredecessors(v);
		for (SeqVertex p : predecessors) {
			
			if (getColor(p) != WHITE) { // shouldnt be white anyway
				if (max_up_depth == null) {
					max_up_depth = p.getNodeDepth();
				}
				else {
					max_up_depth = Math.max(max_up_depth, p.getNodeDepth());
				}
			}
			else {
				//throw new RuntimeException("depth not set yet for node: " + p.getID() + ", pred of " + v.getID());
				if (local_debug) {
					System.err.println("**** " + "depth not set yet for node: " + p.getID() + ", pred of " + v.getID());
				}
			}
		}
		
		
		// redefining depth based on max parent depth value
		Integer node_depth = recommended_depth;
		
		
		if (max_up_depth != null) {
			node_depth = Math.max(max_up_depth + 1, recommended_depth);
		}
		
		
		v.setDepth(node_depth);
		v.setNodeDepth(node_depth);
		
	
		
		setColor(v, DARKGRAY); // depth set based on parents
		
		if (local_debug) {
			System.err.println("*DEPTH(" + v.getID() + ") => " + node_depth + "  descr: " + v.getShortSeqWconnectingIDs(_graph) 
					+ ", recommended depth was: " + recommended_depth);
		}


		*/
		
	
		
		//////////////////////////////////////////
		// Back from visiting children.
		//  Done with all visits.  Exit recursion.
		//////////////////////////////////////////
		
		setColor(v, BLACK);
		_time++;
		// mark the finishing time of this vertex
		_finishing.put(v, _time);
		v.setDFS_FinishingTime(_time);
		
		
		if (local_debug) {
			System.err.println("DFS EXIT for: " + v.getShortSeqWconnectingIDs(_graph));
		}
		
	}
	
	
	/**
	 * Returns the color of the given vertex	
	 */
	public int getColor(SeqVertex v)
	{
		return  _colors.get(v);
	}

	/**
	 * Sets the color of the given vertex	
	 */
	public void setColor(SeqVertex v, Integer color)
	{
		_colors.put(v,color);
	}

	/**
	 * return the finishing time for each vertex
	 * @return
	 */
	public Map<SeqVertex, Number> getFinishing() {
		return _finishing;
	}



}



