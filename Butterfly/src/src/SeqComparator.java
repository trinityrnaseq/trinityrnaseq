import java.io.Serializable;
import java.util.Comparator;


public class SeqComparator implements Comparator<Object>, Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@Override
	public int compare(Object o1, Object o2) {
		int id1 = ((SeqVertex)o1).getID();
		int id2 = ((SeqVertex)o2).getID();
		
		return id1-id2;
	}

}
