#include "graph.h"


// constructor
//currently using default


// destructor
// rlease memory allocated for graph nodes:
Graph::~Graph () {
  for (unsigned int i=0; i < allNodes.size(); i++) {
    delete allNodes[i];
  }
}
  

void Graph::addLinkedNodes (string a, string b) {
  if (VERBOSE >= debug)
    cout << "Adding nodes: " << a << ", " << b << endl;
  
  Graphnode* a_node = getGraphnode(a);
  Graphnode* b_node = getGraphnode(b);
  
  a_node->addLinkedNode(b_node);
  b_node->addLinkedNode(a_node);
  
}


Graphnode* Graph::getGraphnode (string s) {
  if (VERBOSE >= debug)
    cout << "Getting graphnode for : " << s << endl;
  
  map<string, int>::iterator p;
  // see if (s) already exists:
  p = nodeLookup.find(s);
  Graphnode* s_node;
  if (p != nodeLookup.end()) {
    if (VERBOSE >= debug) 
      cout << "Found it." << endl;
    
    int s_index = p->second;
    s_node = allNodes[s_index];
  
  } else {
    if (VERBOSE >= debug)
      cout << "Didn't find it, inserting node." << endl;
  
    s_node = new Graphnode(s);
    int s_index = allNodes.size();
    allNodes.push_back(s_node);
    if (VERBOSE >= debug)
      cout << "Added new node with " << s << " at position: " << s_index << endl;
    nodeLookup[s] = s_index;
  }
  
  return (s_node);
  
}


void Graph::printClusters () {
  vector<string> cluster;
  for (unsigned int i=0; i < allNodes.size(); i++) {
    Graphnode* g = allNodes[i];
    if (g->marked) {
      continue; // already in another cluster
    }
        
    if (!g->marked) {
      
      if (VERBOSE >= debug)
        cout << "Starting graph traversal from node[" << i << "]: " << g->toString() << endl;
      
      traverseGraph(g, cluster);
    }
    // cout << "Cluster : ";
    for (unsigned int j=0; j < cluster.size(); j++) {
      cout << cluster[j] << " ";
    }
    cout << endl;
    cluster.clear(); // init for next cluster
  }
}


void Graph::traverseGraph (Graphnode* g, vector<string>& cluster) {
  if (g->marked) {
    return;
  }
  cluster.push_back(g->getNodename());
  g->marked = true;
  vector<Graphnode*>& linkednodes = g->getLinkedNodes();
  for (unsigned int i=0; i < linkednodes.size(); i++) {
    traverseGraph(linkednodes[i], cluster);
  }
}


 /* Apply the Jaccard Coefficient */
Graph* Graph::applyJaccardCoeff (float coeff) {
  
  Graph* gph = new Graph();
  
  int num_nodes = allNodes.size();
  
  
  for (unsigned int i=0; i < allNodes.size(); i++) {
    
    if (VERBOSE >= info) {
      cerr << "\rJaccard: examining node " << (i+1)  << " of " << num_nodes << "       "; 
    }
    
    Graphnode* g = allNodes[i];
    
    // examine each pair of linked nodes
    vector<Graphnode*>& linkednodes = g->getLinkedNodes();
    for (unsigned int j=0; j < linkednodes.size(); j++) {
      Graphnode* linkednode = linkednodes[j];
      if (calclinkcoeff(g, linkednode) >= coeff) {
        gph->addLinkedNodes(g->getNodename(), linkednode->getNodename());
      }
      
    }
    
  }
  
  return (gph);
  
}


float Graph::calclinkcoeff (Graphnode* a_node, Graphnode* b_node) {
  
  vector<Graphnode*>& a_links = a_node->getLinkedNodes();
  int num_a_links = a_node->numLinkedNodes();
  int num_b_links = b_node->numLinkedNodes();
  
  // determine number shared vertices:
  map<string,bool> b_map = b_node->getLinkedNodeNameMap();
  
  int num_common = 0;
  for (unsigned int i=0; i < a_links.size(); i++) {
    string a_link_name = a_links[i]->getNodename();
    if (b_map.find(a_link_name) != b_map.end()) {
      // found it:
      num_common++;
    }
  }
  
  int total_num_vertices = num_a_links + num_b_links - num_common;
  int total_in_common = num_common + 2; // a and b aren't included in common count
  
  float linkcoeff = (float)total_in_common/total_num_vertices;
  
  if (VERBOSE >= debug) {
    cerr << "Link score between " << a_node->getNodename() << " and " 
         << b_node->getNodename() << " = " << linkcoeff << endl;
  }
  
  return (linkcoeff);
  
  
}
  







