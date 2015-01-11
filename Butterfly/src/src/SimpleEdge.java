
public class SimpleEdge {
	public double _wei;
	public boolean _isInCircle;
	protected int _numberOfLoopsInvolved;
	
	private int from_vertex_id;
	private int to_vertex_id;
	
	protected double _repeat_unroll_weight; // for encouraging unrolled repeat edge use in repeat unrolling
	
	/*
	 * ctor
	 */
	public SimpleEdge (double w, int from_vertex_id, int to_vertex_id)
	{
		_wei = w;
		_isInCircle = false;
		_numberOfLoopsInvolved = 0;
		
		_repeat_unroll_weight = 0;
	
		this.from_vertex_id = from_vertex_id;
		this.to_vertex_id = to_vertex_id;
	
	}

	
	
	public SimpleEdge(SimpleEdge e2, int from_vertex_id, int to_vertex_id) {
		_wei = e2.getWeight();
		_isInCircle = e2.isInCircle();
		_numberOfLoopsInvolved = e2.getNumOfLoopsInvolved();
		this.from_vertex_id = from_vertex_id;
		this.to_vertex_id = to_vertex_id;
	}

	
	
	public double get_repeat_unroll_weight () {
		
		return(_repeat_unroll_weight);
	}

	public void set_repeat_unroll_weight (double w) {
		_repeat_unroll_weight = w;
	}
	
	public double increment_repeat_unroll_weight (double add_w) {
		_repeat_unroll_weight += add_w;
		
		return(_repeat_unroll_weight);
		
	}
	
	/**
	 * @return the number of loops involving this edge
	 */
	public int getNumOfLoopsInvolved() {
		return _numberOfLoopsInvolved;
	}

	public void increaseNumOfLoopsBy1(){
		_numberOfLoopsInvolved++;
	}
	

	public void decreaseNumOfLoopsBy1() {
		_numberOfLoopsInvolved--;		
	}
	/*
	 * return weight
	 */
	public double getWeight()
	{
		return _wei;
	}

	/*
	 * set weight
	 */
	public void setWeight(double w)
	{
		_wei = w;
	}

	
	/**
	 * @return the _isInCircle
	 */
	public boolean isInCircle() {
		return _isInCircle;
	}



	/**
	 * @param isInCircle the _isInCircle to set
	 */
	public void setInCircle(boolean isInCircle) {
		_isInCircle = isInCircle;
	}



	/*
	 * toString
	 */
	public String toString()
	{
		return "Edge(" + this.from_vertex_id + "->" + this.to_vertex_id + ",w:"+_wei+")";
	}





	
//	public boolean equals(Object other)
//	{
//		return other!=null && 
//		_v1.equals(((SimpleEdge) other).getSourceVertex()) &&
//		_v2.equals(((SimpleEdge) other).getTargetVertex());
//	}
//	
//	public int hashCode()
//	{
////		Integer.parseInt("1"+_v1.getID()+"00"+_v2.getID());	
//		return (_v1.getID()+"00"+_v2.getID()).hashCode(); //make a string, and use the string hashcode.
//	}

}
