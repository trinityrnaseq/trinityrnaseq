import java.util.HashMap;
import java.util.Map;

import edu.uci.ics.jung.graph.DirectedSparseGraph;


public class Visitor {
	int _time;
	Map<SeqVertex,Number> _discovery;
	Map<SeqVertex,Number> _finishing;
	Visitor(){
		_time=0;
		_discovery = new HashMap<SeqVertex,Number>();
		_finishing = new HashMap<SeqVertex,Number>();
	}

	public void discoverEdge(SimpleEdge e) {
		// TODO Auto-generated method stub
		
	}

	public void finishEdge(SimpleEdge e) {
		// TODO Auto-generated method stub
		
	}

	public void discoverVertex(SeqVertex v) {
		_discovery.put(v, _time);
		v.setDFS_DiscoveryTime(_time);
//		System.err.println("vertex "+v.getID()+" got discovery time of "+_time);
		_time++;
		
	}

	public void finishVertex(SeqVertex v) {
		_finishing.put(v, _time);
		v.setDFS_FinishingTime(_time);
//		System.err.println("vertex "+v.getID()+" got finishing time of "+_time);
		_time++;
		
	}

	public void discoverGraph(DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		// TODO Auto-generated method stub
		
	}

	public void finishGraph(DirectedSparseGraph<SeqVertex, SimpleEdge> graph) {
		// TODO Auto-generated method stub
		
	}
	
	public Map<SeqVertex,Number> getDiscovery()
	{
		return _discovery;
	}
	
	public Map<SeqVertex,Number> getFinishing()
	{
		return _finishing;
	}

}
