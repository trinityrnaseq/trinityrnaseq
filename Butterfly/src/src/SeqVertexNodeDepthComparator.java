import java.util.Comparator;

public class SeqVertexNodeDepthComparator implements Comparator {

	@Override
	public int compare(Object o1, Object o2) {
		int f1 = ((SeqVertex)o1).getDepth();
		int f2 = ((SeqVertex)o2).getDepth();
		if( f1 < f2 )
			return -1;
		else if( f1 > f2 )
			return 1;
		else
			return 0;
	}
	
	
}

