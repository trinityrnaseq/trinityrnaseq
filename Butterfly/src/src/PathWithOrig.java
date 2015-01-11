import java.util.*;


public class PathWithOrig {

	private String pathNodeID = null;
	private List<Integer> vertex_id_list;
	private List<Integer> orig_vertex_id_list;
	
	public PathWithOrig(List<Integer> path) {
		
		vertex_id_list = new ArrayList<Integer>();
		vertex_id_list.addAll(path);
		
		orig_vertex_id_list = new ArrayList<Integer>();
		for (Integer node_id : vertex_id_list) {
			Integer orig_id = SeqVertex.retrieveSeqVertexByID(node_id).getOrigButterflyID();
			orig_vertex_id_list.add(orig_id);
		}
		
	}
	
	public PathWithOrig(List<Integer> vertex_id_list, List<Integer> orig_vertex_id_list) {
		this.vertex_id_list = vertex_id_list;
		this.orig_vertex_id_list = orig_vertex_id_list;
	}
	
	
	public PathWithOrig(String pathNodeID, List<Integer> vertex_id_list,
			List<Integer> orig_vertex_id_list) {
		
		this(vertex_id_list, orig_vertex_id_list);
		
		this.pathNodeID = pathNodeID;
		
	}

	public List<Integer> getVertexList() {
		
		return(vertex_id_list);
	}
	
	public List<Integer> getOrigVertexList() {
		return(orig_vertex_id_list);
	}
	
	
	public int size () {
		return(vertex_id_list.size());
	}

	public PathWithOrig align_path_by_orig_id(PathWithOrig template_pwo) {
		
		//System.err.println("\nComparing read: " + this + " to Template: " + template_pwo);
		
		List<Integer> template_orig_vertex_list = template_pwo.getOrigVertexList(); 
		
		for (int i = 0; i < template_orig_vertex_list.size(); i++) {
			if (this.orig_vertex_id_list.get(0).equals(template_orig_vertex_list.get(i))) {
				
				List<Integer> restructured_path = new ArrayList<Integer>();
				int matches = 0;
				
				// have a first position match.  See if it extends to a complete match:
				for (int this_path_pos=0, template_path_pos = i; 
						template_path_pos < template_orig_vertex_list.size() 
						 && this_path_pos < this.orig_vertex_id_list.size();
						template_path_pos++, this_path_pos++) {
				
					Integer template_orig_id = template_orig_vertex_list.get(template_path_pos);
					Integer this_orig_id = this.orig_vertex_id_list.get(this_path_pos);
					
					
					if (this_orig_id.equals(template_orig_id)) {
						matches++;
						restructured_path.add(template_pwo.vertex_id_list.get(template_path_pos));
						//System.err.println("" + this_path_pos + "=" + this_orig_id
						//		 + " vs. " + template_path_pos + "=" + template_orig_id + " matches.");
					}
					else {
						restructured_path.add(template_pwo.vertex_id_list.get(template_path_pos));
						//System.err.println("" + this_path_pos + "=" + this_orig_id
						//		 + " vs. " + template_path_pos + "=" + template_orig_id + " MISmatches.");
						break;
					}
				}
				
				//System.err.println("Got match count: " + matches + " vs. path size " + this.orig_vertex_id_list.size());
				
				if (matches == this.orig_vertex_id_list.size()) {
					// got it.
					//System.err.println("** restructured path OK");
					return(new PathWithOrig(restructured_path, this.orig_vertex_id_list));
					
				}
				
			}
		}
		
		//System.err.println("*!* no proper restrucuring of path.");		
		return(null); // not found.
		
	}
	
	public String toString() {
		String text = "";
		if (this.pathNodeID != null) {
			text += this.pathNodeID + " ";
		}
		
		text += "CurrVert:" + this.vertex_id_list + ", OrigVert:" + this.orig_vertex_id_list;
		
		return(text);
	}

	public void update_vertex_list(List<Integer> updated_path) {
		this.vertex_id_list = updated_path;
	}
	
}
