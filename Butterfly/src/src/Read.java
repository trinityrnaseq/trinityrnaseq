import java.util.List;



public class Read {
	private String _name;
	private String _seq;
	private LocInGraph _startingV;
	private Integer _startInRead;
	private Integer _endInRead;
	private List<Integer> _path_ids;
	/**
	 * Ctor
	 * @param name
	 * @param seq
	 * @param endInRead 
	 * @param startInRead 
	 * @param toV 
	 */
	public Read(String name, String seq, LocInGraph startingV, Integer startInRead, Integer endInRead, List<Integer> pathIds) {
		super();
		init( name,  seq,  startingV,  startInRead,  endInRead,  pathIds);
	}

	/** 
	 * init function, to use after the default ctor.
	 * @param name
	 * @param seq
	 * @param startingV
	 * @param startInRead
	 * @param endInRead
	 * @param pathIds
	 */
	public void init (String name, String seq, LocInGraph startingV, 
			Integer startInRead, Integer endInRead, List<Integer> pathIds) {	
		_name = name;
		_seq = seq;
		_startingV = startingV;
		_startInRead = startInRead;
		_endInRead = endInRead;
		_path_ids = pathIds;
	}

	/** 
	 * empty ctor
	 */
	public Read() {
		super();
		_name = "";
		_seq = "";
		_startingV = new LocInGraph(-1, -1);
		_startInRead = -1;
		_endInRead = -1;
		_path_ids = null;
	}


	/**
	 * Return the name of this read
	 * @return
	 */
	public String getName() {
		return _name;
	}

	/**
	 * Return the seq of this read
	 * @return
	 */
	public String getSeq() {
		return _seq;
	}

	/**
	 * Return the starting id of this read
	 * @return
	 */
	public LocInGraph getStartingV() {
		return _startingV;
	}

	

	/**
	 * Return the first index of the match of the read to the graph
	 * @return
	 */
	public Integer getStartInRead() {
		return _startInRead;
	}

	/**
	 * Return the last index of the match of the read to the graph
	 * @return
	 */
	public Integer getEndInRead() {
		return _endInRead;
	}

	/**
	 * Return the list of all ids of this read
	 * @return
	 */
	public List<Integer> getPathIDs() {
		return _path_ids;
	}

	/**
	 * toString
	 */
	@Override
	public String toString() {
		return "Read [_name=" + _name + ", _path_ids=" + _path_ids + ", _seq="
				+ _seq + "]";
	}
	/**
	 * hashcode
	 */
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((_name == null) ? 0 : _name.hashCode());
		result = prime * result
				+ ((_path_ids == null) ? 0 : _path_ids.hashCode());
		result = prime * result + ((_seq == null) ? 0 : _seq.hashCode());
		result = prime * result
				+ ((_startingV == null) ? 0 : _startingV.hashCode());
		return result;
	}
	
	/**
	 * equals
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Read other = (Read) obj;
		if (_name == null) {
			if (other._name != null)
				return false;
		} else if (!_name.equals(other._name))
			return false;
		if (_path_ids == null) {
			if (other._path_ids != null)
				return false;
		} else if (!_path_ids.equals(other._path_ids))
			return false;
		if (_seq == null) {
			if (other._seq != null)
				return false;
		} else if (!_seq.equals(other._seq))
			return false;
		if (_startingV == null) {
			if (other._startingV != null)
				return false;
		} else if (!_startingV.equals(other._startingV))
			return false;
		return true;
	}
	
	

}
