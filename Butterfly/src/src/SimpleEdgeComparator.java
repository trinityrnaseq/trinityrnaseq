import java.util.Comparator;


public class SimpleEdgeComparator implements Comparator<Object> {

	@Override
	public int compare(Object o1, Object o2) {
		double f1 = ((SimpleEdge)o1).getWeight();
		double f2 = ((SimpleEdge)o2).getWeight();
		if( f1 < f2 )
			return 1;
		else if( f1 > f2 )
			return -1;
		else
			return 0;

	
	}

}
