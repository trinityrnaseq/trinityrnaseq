
public class SimplePathNodeEdge {

	float weight;
	String from_PathNodeID;
	String to_PathNodeID;
	
	int _num_loops;
	
	
	public SimplePathNodeEdge(float weight, String from_pathNodeID, String to_pathNodeID) {
		
		this.weight = weight;
		this.from_PathNodeID = from_pathNodeID;
		this.to_PathNodeID = to_pathNodeID;
		
		this._num_loops = 0;
		
	}
	
	
	public String toString () {
		
		return("PNEdge(" + from_PathNodeID + "->" + to_PathNodeID + ")");
	}


	public void increaseNumOfLoopsBy1() {
		
		_num_loops++;
		
	}


	public int getNumOfLoopsInvolved() {
		return(_num_loops);
	}


	public void decreaseNumOfLoopsBy1() {
		_num_loops--;
		
	}
	
	
}
