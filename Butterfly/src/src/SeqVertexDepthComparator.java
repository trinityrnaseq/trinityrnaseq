import java.util.Comparator;


public class SeqVertexDepthComparator implements Comparator<Object> {

	@Override
	public int compare(Object o1, Object o2) {
		int d1 = ((SeqVertex)o1).getDepth();
		int d2 = ((SeqVertex)o2).getDepth();
		if( d1 > d2 )
			return 1;
		else if( d1 < d2 )
			return -1;
		else
			return 0;
	}

}
