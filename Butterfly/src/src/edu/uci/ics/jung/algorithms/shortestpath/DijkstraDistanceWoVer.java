package edu.uci.ics.jung.algorithms.shortestpath;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.graph.Graph;


public class DijkstraDistanceWoVer<V, E> extends DijkstraDistance<V, E> {

	public DijkstraDistanceWoVer(Graph<V, E> g) {
//		super(g,true);
		super(g,false);

	}


	public Number getDistanceWoVer(V source, V target,V verToExclude) 
	{
		if (g.containsVertex(target) == false)
			throw new IllegalArgumentException("Specified target vertex " + 
					target + " is not part of graph " + g);
		if (g.containsVertex(source) == false)
			throw new IllegalArgumentException("Specified source vertex " + 
					source + " is not part of graph " + g);

		Set<V> targets = new HashSet<V>();
		targets.add(target);
		Map<V,Number> distanceMap = getDistanceMap(source, targets,verToExclude);
		return distanceMap.get(target);
	}


	public Map<V,Number> getDistanceMap(V source, Collection<V> targets, V verToExclude)
	{
		if (g.containsVertex(source) == false)
			throw new IllegalArgumentException("Specified source vertex " + 
					source + " is not part of graph " + g);
		if (targets.size() > max_targets)
			throw new IllegalArgumentException("size of target set exceeds maximum " +
					"number of targets allowed: " + this.max_targets);

		Map<V,Number> distanceMap = 
			singleSourceShortestPath(source, targets, 
					Math.min(g.getVertexCount(), max_targets),verToExclude);
		if (!cached)
			reset(source);

		return distanceMap;
	}


	protected LinkedHashMap<V,Number> singleSourceShortestPath(V source, Collection<V> targets, int numDests, V verToExclude)
	{
		SourceData sd = getSourceData(source);

		Set<V> to_get = new HashSet<V>();
		if (targets != null) {
			to_get.addAll(targets);
			Set<V> existing_dists = sd.distances.keySet();
			for(V o : targets) {
				if (existing_dists.contains(o))
					to_get.remove(o);
			}
		}

		// if we've exceeded the max distance or max # of distances we're willing to calculate, or
		// if we already have all the distances we need, 
		// terminate
		if (sd.reached_max ||
				(targets != null && to_get.isEmpty()) ||
				(sd.distances.size() >= numDests))
		{
			return sd.distances;
		}

		while (!sd.unknownVertices.isEmpty() && (sd.distances.size() < numDests || !to_get.isEmpty()))
		{
			Map.Entry<V,Number> p = sd.getNextVertex();
			V v = p.getKey();
			double v_dist = ((Double)p.getValue()).doubleValue();
			sd.dist_reached = v_dist;
			to_get.remove(v);
			if ((sd.dist_reached >= this.max_distance) || (sd.distances.size() >= max_targets))
			{
				sd.reached_max = true;
				break;
			}

			for (E e : getEdgesToCheck(v) )
			{
				for (V w : g.getIncidentVertices(e))
				{
					if (!w.equals(verToExclude) && !sd.distances.containsKey(w))
					{
						double edge_weight = nev.transform(e).doubleValue();
						if (edge_weight < 0)
							throw new IllegalArgumentException("Edges weights must be non-negative");
						double new_dist = v_dist + edge_weight;
						if (!sd.estimatedDistances.containsKey(w))
						{
							sd.createRecord(w, e, new_dist);
						}
						else
						{
							double w_dist = ((Double)sd.estimatedDistances.get(w)).doubleValue();
							if (new_dist < w_dist) // update tentative distance & path for w
								sd.update(w, e, new_dist);
						}
					}
				}
			}
			//            // if we have calculated the distance to the target, stop
			//            if (v == target)
			//                break;

		}
		return sd.distances;
	}
}

