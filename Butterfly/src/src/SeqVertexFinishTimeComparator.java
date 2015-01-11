import java.util.Comparator;


public class SeqVertexFinishTimeComparator implements Comparator<Object> {
	@Override
	public int compare(Object o1, Object o2) {
		int f1 = ((SeqVertex)o1).getDFS_FinishingTime();
		int f2 = ((SeqVertex)o2).getDFS_FinishingTime();
		if( f1 < f2 )
			return 1;
		else if( f1 > f2 )
			return -1;
		else
			return 0;
	}


}

