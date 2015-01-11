import java.util.List;


public class ScoredPath {


	List<PairPath> paths;
	int score = 0;
	boolean path_extended = false;

	public ScoredPath (List<PairPath> paths, int score) {
		this.paths = paths;
		this.score = score;
	}
	
	
	public String toString() {
		String ret = "Score: " + score + " PPaths: " + paths;
		return(ret);
	}
	
}



