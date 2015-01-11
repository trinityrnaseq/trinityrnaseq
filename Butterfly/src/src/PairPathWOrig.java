import java.util.*;


public class PairPathWOrig {

	private PathWithOrig path1;
	private PathWithOrig path2;
	
	public PairPathWOrig(PairPath pp) {
		
		path1 = new PathWithOrig(pp.getPath1());
		
		path2 = null;
		if (pp.hasSecondPath()) {
			path2 = new PathWithOrig(pp.getPath2());
		}
		
		
	}
	
	public PairPathWOrig(PathWithOrig p1, PathWithOrig p2) {
		path1 = p1;
		path2 = p2;
	}
	
	public PairPathWOrig(PathWithOrig p1) {
		path1 = p1;
		path2 = null;
	}
	
	public PathWithOrig getPath1() {
		return(path1);
	}
	
	public PathWithOrig getPath2() {
		return(path2);
	}
	
	public boolean hasPath2() {
		return(path2 != null);
	}
	
	public PairPath getPairPath() {
		
		PairPath pp;
		
		if (path2 == null) {
			pp = new PairPath(path1.getVertexList());
		}
		else {
			pp = new PairPath(path1.getVertexList(), path2.getVertexList());
		}
		
		return(pp);
	}
	
	public int size() {
		int s = path1.size();
		if (path2 != null) {
			s = Math.max(s,path2.size());
		}
		return(s);
	}

	public PairPathWOrig restructure_according_to_repeat_path(PathWithOrig template_pwo) {
		
		PathWithOrig pwo_restructured = this.path1.align_path_by_orig_id(template_pwo);
		if (pwo_restructured == null) {
			return(null);
		}
		if (! this.hasPath2()) {
			return (new PairPathWOrig(pwo_restructured));
		}
		
		// must have path2
		PathWithOrig pwo2_restructured = this.path2.align_path_by_orig_id(template_pwo);
		if (pwo2_restructured == null) {
			return(null);
		}
		else {
			return(new PairPathWOrig(pwo_restructured, pwo2_restructured));
		}
		
	}
	
	public boolean equals (PairPathWOrig ppwo) {
		if (! this.path1.getVertexList().equals(ppwo.path1.getVertexList())) {
			return(false);
		}
		if (this.hasPath2() != ppwo.hasPath2()) {
			return(false);
			
		}
		if (this.hasPath2() && ! this.path2.getVertexList().equals(ppwo.path2.getVertexList()) ) {
			return(false);
		}
		
		// must be the same
		return(true);
		
	}
	
	public String toString() {
		
		// path 1
		String text = "p1: {" +  path1.toString() + "}";
		if (this.hasPath2()) {
			text += ", p2: {" + path2.toString() + "}";
		}
		
		return(text);
	}
	
	
	
}
