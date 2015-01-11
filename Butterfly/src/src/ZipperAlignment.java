import java.util.HashMap;


public class ZipperAlignment {
	
	//private int kmer_size = 8;
	
	public static AlignmentStats doZipperAlignment(String nameA, String seqA, String nameB, String seqB) {
		
		AlignmentStats anchorLeftStats = doZipperAlignmentAnchorLeft(nameA, seqA, nameB, seqB);
		
		AlignmentStats anchorRightStats = doZipperAlignmentAnchorRight(nameA, seqA, nameB, seqB);
		
		
		if (anchorLeftStats.matches > anchorRightStats.matches) {
			return(anchorLeftStats);
		}
		else {
			return(anchorRightStats);
		}
		
	}
		
	
	public static AlignmentStats doZipperAlignmentAnchorLeft(String nameA, String seqA, String nameB, String seqB) {
		
		int shorter_len, longer_len;
		if (seqA.length() < seqB.length()) {
			shorter_len = seqA.length();
			longer_len = seqB.length();
		}
		else {
			shorter_len = seqB.length();
			longer_len = seqA.length();
		}
		
		AlignmentStats as = new AlignmentStats();
		as.alignedBaseCounter = new HashMap<String,Integer>();
	    as.alignedBaseCounter.put(nameA, new Integer(shorter_len));
	    as.alignedBaseCounter.put(nameB, new Integer(shorter_len));
		
		
		for (int i = 0; i < shorter_len; i++) {
			if (seqA.charAt(i) == seqB.charAt(i)) {
				as.matches++;
			}
			else {
				as.mismatches++;
			}
		}
		as.right_gap_length = longer_len - shorter_len;
		as.total_not_matched = as.right_gap_length;		
		
		return(as);
	}
	
public static AlignmentStats doZipperAlignmentAnchorRight(String nameA, String seqA, String nameB, String seqB) {
		
		int shorter_len, longer_len;
		if (seqA.length() < seqB.length()) {
			shorter_len = seqA.length();
			longer_len = seqB.length();
		}
		else {
			shorter_len = seqB.length();
			longer_len = seqA.length();
		}
	
		AlignmentStats as = new AlignmentStats();
		as.alignedBaseCounter = new HashMap<String,Integer>();
	    as.alignedBaseCounter.put(nameA, new Integer(shorter_len));
	    as.alignedBaseCounter.put(nameB, new Integer(shorter_len));
		
		
		for (int i = shorter_len-1, j = longer_len-1; i >= 0 && j >= 0; i--, j--) {
			if (seqA.charAt(i) == seqB.charAt(i)) {
				as.matches++;
			}
			else {
				as.mismatches++;
			}
		}
		as.left_gap_length = longer_len - shorter_len;
		as.total_not_matched = as.left_gap_length;		
		
		return(as);
	}
	
	
	

}
