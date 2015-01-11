#ifndef __GRAPHNODE__
#define __GRAPHNODE__

#include <vector>
#include <string>
#include <map>
#include "common.h"

using namespace std;

class Graphnode {
	public:
	
	string nodeName;
	vector<Graphnode*> linkedNodes;
	bool marked;
	
    // constructor
	Graphnode (string s);
	
    
    // get name of node
	string getNodename ();
	

    // add a linked node
	void addLinkedNode (Graphnode* a);

	
    // get number of linked nodes.
	int numLinkedNodes ();
    	
    
    // describe node
    string toString ();
    	
    
    // get list of linked nodes by their names
	map<string,bool> getLinkedNodeNameMap ();
    
    // get the actual linked nodes list
	vector<Graphnode*>& getLinkedNodes ();
};

#endif
