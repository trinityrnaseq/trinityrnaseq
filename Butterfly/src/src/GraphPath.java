import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * a class representing a path in the graph
 * @author moran
 *
 */
public class GraphPath {

	protected List<SeqVertex> _path;
	protected double _totalSupport;
	protected int _size;
	protected boolean _isCircular; 
	
	
	public GraphPath(List<SeqVertex> path) {
		_path = path;
		_size = 0;
		_totalSupport = 0;
		_isCircular = false;
		for (SeqVertex v : path)
		{
			addToSupport(v.getWeightSum());
			addToSize(v.getSize());
		}
	}

	public GraphPath() {
		_path = new ArrayList<SeqVertex>();
		_totalSupport = 0;
		_size = 0;
		_isCircular = false;
	}


	/**
	 * @return
	 * @see java.util.List#iterator()
	 */
	public Iterator<SeqVertex> iterator() {
		return _path.iterator();
	}

	@Override
	public String toString() {
		return "GraphPath [" + _path + "]("+_totalSupport+")";
			
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((_path == null) ? 0 : _path.hashCode());
		long temp;
		temp = Double.doubleToLongBits(_totalSupport);
		result = prime * result + (int) (temp ^ (temp >>> 32));
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
		GraphPath other = (GraphPath) obj;
		if (_path == null) {
			if (other._path != null)
				return false;
		} else if (!_path.equals(other._path))
			return false;
		if (Double.doubleToLongBits(_totalSupport) != Double
				.doubleToLongBits(other._totalSupport))
			return false;
		return true;
	}

	/**
	 * @return the _path
	 */
	public List<SeqVertex> getPath() {
		return _path;
	}

	public void addToPath(SeqVertex v)
	{
		_path.add(v);
	}
	
	/**
	 * @return the _support
	 */
	public double getSupport() {
		return _totalSupport/_size;
	}

	/**
	 * @param support the _support to set
	 */
	public void setSupport(double support) {
		_totalSupport = support;
	}

	
	public void addToSupport(double weightSum) {
		_totalSupport += weightSum;		
	}

	public void addToSize(int size) {
		_size += size;
	}

	/**
	 * @return the _isCircular
	 */
	public boolean isCircular() {
		return _isCircular;
	}

	/**
	 * set this path as circular
	 */
	public void setCircular() {
		_isCircular = true;
	}



}
