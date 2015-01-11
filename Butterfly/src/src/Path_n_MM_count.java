import java.util.ArrayList;
import java.util.List;
import java.util.Vector;


public class Path_n_MM_count {

	Integer mismatch_count = 0;
	Vector<Integer> path;
	Vector<Vertex_path_position> path_read_positions;
	
	

	public Path_n_MM_count() {
		path = new Vector<Integer>();
		path_read_positions = new Vector<Vertex_path_position>();
		
	}

	public Path_n_MM_count(Path_n_MM_count orig) {
		
		this();
		this.mismatch_count = orig.mismatch_count;
		this.path.addAll(orig.path);
		this.path_read_positions.addAll(orig.path_read_positions);
		
	}
	
	
	public Path_n_MM_count(Integer node, Integer mm, Integer end_vertex_base, Integer end_read_base) {

		path = new Vector<Integer>();
		path_read_positions = new Vector<Vertex_path_position>();
		
		mismatch_count += mm;
		path.add(0, node);
	
		Vertex_path_position vpp = new Vertex_path_position(node, end_vertex_base, end_read_base);
		path_read_positions.add(0, vpp);
	}

	public void add_path_n_mm (Integer node, Integer mm, Integer end_vertex_base, Integer end_read_base) {
		path.add(0, node);
		mismatch_count += mm;
		
		Vertex_path_position vpp = new Vertex_path_position(node, end_vertex_base, end_read_base);
		path_read_positions.add(0, vpp);
	}
	
	
	public String toString() {
		String text = "PATH_N_MM_COUNT: mismatches=" + mismatch_count + ", path= " + path + " positions:  ";
		
		for (Vertex_path_position vpp: path_read_positions) {
			text += "[" + vpp.vertex_id + ",nodeEnd:" + vpp.vertex_base_end + ",readEnd:" + vpp.read_base_end + "] ";
		}
		
		return(text);
		
	}
	
	public List<Integer> get_trimmed_path (Integer min_end_length) {
		 
		Vector<Vertex_path_position> vpps = new Vector<Vertex_path_position>(path_read_positions);
		
		while ( (! vpps.isEmpty()) && vpps.elementAt(0).read_base_end < min_end_length) {
			
			vpps.removeElementAt(0);
		}
		
		while ( (! vpps.isEmpty()) && vpps.elementAt(vpps.size()-1).vertex_base_end < min_end_length) {
			vpps.removeElementAt(vpps.size()-1);
		}
		
		List<Integer> trimmed_path = new ArrayList<Integer>();
		
		for (Vertex_path_position vpp : vpps) {
			trimmed_path.add(vpp.vertex_id);
		}
		
		return(trimmed_path);
	}
	
	
	
	class Vertex_path_position {
		
		Integer vertex_id;
		Integer vertex_base_end;
		Integer read_base_end;
		
		public Vertex_path_position(Integer vertex_id, Integer vertex_base_end, Integer read_base_end) {
			
			this.vertex_id = vertex_id;
			this.vertex_base_end = vertex_base_end;
			this.read_base_end = read_base_end;
		}
	}
	
	

}
