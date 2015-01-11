import java.util.*;


public class PasaVertex {

	final PairPath pp;
	int readSupport = 0;
	
	int vertex_score = 0;
	
	public static int max_top_paths_to_store = 1;
	
	List<ScoredPath> fromPaths;
	List<ScoredPath> toPaths;
	
	int num_contained = 0;
	
	private Comparator<ScoredPath> sp_comparator = new Comparator<ScoredPath>() {
		
		public int compare (ScoredPath spA, ScoredPath spB) {
			
			// want it to sort descendingly on score
			
			if (spA.score < spB.score) {
				return(1);
			}
			else if (spA.score > spB.score) {
				return(-1);
			}
			else {
				return(0);
			}
		}
		
	};
	
	
	public PasaVertex (final PairPath p, int readSupport) {
		
		pp = p;
		this.readSupport = readSupport;	
		
		
		fromPaths = new ArrayList<ScoredPath>();
		toPaths = new ArrayList<ScoredPath>();
		
		List<PairPath> path = new ArrayList<PairPath>();
		path.add(p);
		
		ScoredPath sp = new ScoredPath(path, num_contained + 1);
		
		fromPaths.add(sp);
		toPaths.add(sp);
	}
	
	
	public final List<ScoredPath> get_toPaths () {
		return(toPaths);
	}
	
	public final List<ScoredPath> get_fromPaths() {
		return(fromPaths);
	}
	
	
	public void push_toPaths(ScoredPath sp) {
		
		push_path_list(this.toPaths, sp);
		
		
	}
	
	public void push_fromPaths(ScoredPath sp) {
		
		push_path_list(this.fromPaths, sp);
		
	}
	
	
	private void push_path_list (List<ScoredPath> sp_list, ScoredPath sp) {
	
		//System.out.println("SP_LIST size before: " + sp_list.size());
		
		sp_list.add(sp);
		
		if (sp_list.size() > PasaVertex.max_top_paths_to_store) {
			Collections.sort(sp_list, this.sp_comparator);
		
			sp_list.retainAll(sp_list.subList(0, max_top_paths_to_store));
		}
		
		//System.out.println("SP_LIST size after: " + sp_list.size());
	
	}
	
	
	public final ScoredPath get_highest_scoring_toPath () {
		
		ScoredPath best = null;
		for (ScoredPath sp : this.get_toPaths()) {
			if (best == null || best.score < sp.score)
				best = sp;
		}
		return(best);
	}
	
	
public final List<ScoredPath> get_all_highest_scoring_toPath () {
		
		List<ScoredPath> best_paths = new ArrayList<ScoredPath>();
		
		for (ScoredPath sp : this.get_toPaths()) {
			if (best_paths.isEmpty()) {
				best_paths.add(sp);
			}
			else if (sp.score == best_paths.get(0).score) {
				best_paths.add(sp);
			}
			else if (sp.score > best_paths.get(0).score) {
				best_paths.clear();
				best_paths.add(sp);
			}
			
		}
		return(best_paths);
	}
	
	
	
	public final ScoredPath get_highest_scoring_fromPath () {
		ScoredPath best = null;
		for (ScoredPath sp : this.get_fromPaths()) {
			if (best == null || best.score < sp.score)
				best = sp;
		}
		return(best);
	}
	
	public final List<ScoredPath> get_all_highest_scoring_fromPath () {

		List<ScoredPath> best_paths = new ArrayList<ScoredPath>();

		for (ScoredPath sp : this.get_fromPaths()) {
			if (best_paths.isEmpty()) {
				best_paths.add(sp);
			}
			else if (sp.score == best_paths.get(0).score) {
				best_paths.add(sp);
			}
			else if (sp.score > best_paths.get(0).score) {
				best_paths.clear();
				best_paths.add(sp);
			}

		}
		return(best_paths);
	}
	
	
	
	public String toString () {
		String ret = "PasaVertex.  PairPath: " + this.pp + "\n";
		ret += "From ScorePaths:\n";
		
		for (ScoredPath sp : this.get_fromPaths()) {
			ret += "\tScore: " + sp.score + ", pp: " + sp.paths + "\n";
			
		}
		
		ret += "To ScorePaths:\n";
		for (ScoredPath sp : this.get_toPaths()) {
			ret += "\tScore: " + sp.score + ", pp: " + sp.paths + "\n\n";
		}
		
		return(ret);
	}
	
	
}
