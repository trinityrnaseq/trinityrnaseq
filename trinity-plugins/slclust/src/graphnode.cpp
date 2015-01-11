#include "graphnode.h"


// constructor
Graphnode::Graphnode (string s) {
  this->nodeName = s;
  marked = false;
}



string Graphnode::getNodename () {
  return (this->nodeName);
}

void Graphnode::addLinkedNode (Graphnode* a) {
  bool found = false;
  string a_name = a->getNodename();
  for (unsigned int i=0; i < linkedNodes.size(); i++) {
    string s = linkedNodes[i]->getNodename();
    if (s.compare(a_name) == 0) {
      found = true;
      break;
    }
  }
  if (! found) {
    linkedNodes.push_back(a);
  }
}


	
int Graphnode::numLinkedNodes () {
  return (linkedNodes.size());
}


string Graphnode::toString () {
  string marked = (this->marked) ? "true" : "false";
  string ret = "(" + nodeName + ", marked: " + marked + ") linked to the following: ";
  for (unsigned int i=0; i < linkedNodes.size(); i++) {
    Graphnode* g = linkedNodes[i];
    string othermarked = (g->marked) ? "true" : "false";
    ret += "(" + g->getNodename() + ", marked: " + othermarked + ") ";
    if (i < linkedNodes.size() - 1) {
      ret += ", ";
    }
  }
  
  return (ret);
}



map<string,bool> Graphnode::getLinkedNodeNameMap () {
  map<string,bool> m;
  for (unsigned int i=0; i < linkedNodes.size(); i++) {
    m [ linkedNodes[i]->getNodename() ] = true;
  }
  return (m);
}

vector<Graphnode*>& Graphnode::getLinkedNodes () {
  return (linkedNodes);
}
