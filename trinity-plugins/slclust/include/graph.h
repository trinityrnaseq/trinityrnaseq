#ifndef __GRAPH__
#define __GRAPH__

#include <vector>
#include <map>
#include <string>
#include "graphnode.h"
#include <iostream>
#include "common.h"

using namespace std;

class Graph {
  
 public:
  
  // all vertices
  vector<Graphnode*> allNodes;
  
  // map of vertex name to index in above vector
  map<string,int> nodeLookup;
  
  // destructor
  ~Graph(); // release memory allocated for graphnodes.
  
  
  // add a pair of linked vertices
  void addLinkedNodes(string a, string b);
  
  // helper function, get a graph node.
  Graphnode* getGraphnode (string s);
  
  // generate cluster output
  void printClusters ();
  
  
  // populates cluster list by finding connected components in graph
  void traverseGraph (Graphnode* g, vector<string>& cluster);
    
  
  /* Apply the Jaccard Coefficient */
  Graph* applyJaccardCoeff (float coeff);

  /* compare two nodes to get their link score (jaccard coeff) */
  float calclinkcoeff (Graphnode* a_node, Graphnode* b_node);
  
};	


#endif		
