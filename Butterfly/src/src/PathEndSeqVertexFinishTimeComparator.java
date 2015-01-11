import java.util.Comparator;


public class PathEndSeqVertexFinishTimeComparator implements Comparator<Object> {
	public int compare(Object o1, Object o2) {
		GraphPath p1 = (GraphPath) o1;
		GraphPath p2 = (GraphPath) o2;

		SeqVertex v1 = p1.getPath().get(p1.getPath().size()-1);
		SeqVertex v2 = p2.getPath().get(p2.getPath().size()-1);
		
		int f1 = v1.getDFS_FinishingTime();
		int f2 = v2.getDFS_FinishingTime();
		if( f1 < f2 )
			return 1;
		else if( f1 > f2 )
			return -1;
		else
		{
			double supp1 = p1.getSupport();
			double supp2 = p2.getSupport();
			if (supp1 < supp2)
				return -1;
			else if (supp1 > supp2)
				return 1;
			else
				return 0;
		}
	}


}

