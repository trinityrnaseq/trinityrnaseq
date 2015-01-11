import java.util.Arrays;


public class BipartiteMatching {
	private int leftSize;
	private int rightSize;
	private int[][] graph;
	private int[] matchingLeft; //left vertex index i is matched to
	private int[] matchingRight; //right vertex index i is matched to
	private boolean[] visited;
	
	public BipartiteMatching(int left, int right, int[][] someGraph)
	{
		leftSize = left;
		rightSize = right;
		graph = someGraph;
		visited = new boolean[graph.length];
		matchingLeft = new int[rightSize];
		matchingRight = new int[leftSize];
	}
	public void printRightMatching()
	{
		for(int i = 0; i < leftSize; i++)
		{
			System.out.println("Match node " + i + " to node " + matchingRight[i]);
		}
	}
	public int[] getLeftMatching()
	{
		return matchingLeft;
	}
	public int[] getRightMatching()
	{
		return matchingRight;
	}
	public void maxMatching()
	{
		 //index is the vertex on the left side and value is its matching vertex on right side
		Arrays.fill(matchingLeft, -1);
		Arrays.fill(matchingRight, -1);
		for(int i = 0; i < leftSize; i++)
		{
			Arrays.fill(visited, false);
			if(bipartiteMatch(i))
				continue;
		}
	}
	private boolean bipartiteMatch(int v)
	{
		for(int j = 0; j < rightSize; j++)
		{
			if(graph[v][j] == 0 || visited[j]) //if the edge doesn't exist or you've already visited this node
				continue;
			visited[j] = true;
			if(matchingRight[j] == -1 || bipartiteMatch(matchingRight[j]))
			{
				matchingLeft[v] = j;
				matchingRight[j] = v;
				//matching[j] = v;
				return true;
			}
		}
		return false;
	}
}
