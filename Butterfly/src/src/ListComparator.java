import java.io.Serializable;
import java.util.List;
import java.util.Comparator;
import java.util.ListIterator;


public class ListComparator implements Comparator<Object>, Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	@SuppressWarnings("unchecked")
	@Override
	public int compare(Object o1, Object o2) {

		int o1_s = ((List)o1).size();	
		int o2_s = ((List)o2).size();
		
//		System.err.println("comparing "+((List)o1) +" to "+((List)o2));
		
		if (o1_s < o2_s)
			return -1;
		else if (o1_s > o2_s)
			return 1;
		else 
		{
			ListIterator i1 = ((List)o1).listIterator();
			ListIterator i2 = ((List)o2).listIterator();
			SeqVertex v1=null, v2=null;
			while (i1.hasNext() && i2.hasNext())
			{
				v1 = (SeqVertex)i1.next();
				v2 = (SeqVertex)i2.next();
				
				if (!v1.equals(v2))
					return v1.getID()-v2.getID();
			}
			
		}
		return 0;
	}
}
