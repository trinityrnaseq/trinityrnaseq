/**
 * This class hold a location in the graph, which is comprised of the node ID, and the location within the node
 * @author moran
 *
 */
public class LocInGraph {
	private int _nodeID;
	private int _indexInNode;
	
	
	public LocInGraph(int nodeID, int indexInNode) {
		super();
		_nodeID = nodeID;
		_indexInNode = indexInNode;
	}


	public int getNodeID() {
		return _nodeID;
	}


	public int getIndexInNode() {
		return _indexInNode;
	}


	public boolean sameNode(LocInGraph other) {
		return getNodeID()==other.getNodeID();
	}


	@Override
	public String toString() {
		return "LocInGraph [_indexInNode=" + _indexInNode + ", _nodeID="
				+ _nodeID + "]";
	}


	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + _indexInNode;
		result = prime * result + _nodeID;
		return result;
	}


	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		LocInGraph other = (LocInGraph) obj;
		if (_indexInNode != other._indexInNode)
			return false;
		if (_nodeID != other._nodeID)
			return false;
		return true;
	}
	
	
}
