#include "KtreeNode.hpp"
#include "stacktrace.hpp"
#include <iostream>
#include <sstream>


KtreeNode::KtreeNode(char nodeChar, long count) {

    this->nodeChar = nodeChar;
    this->count = count;

    for (int i = 0; i < 4; i++) {
        children[i] = 0;
    }
}



bool KtreeNode::has_child (char nodeChar) {
    
    int pos = char_to_position(nodeChar);
    
    if (children[pos] > 0) {
        return(true);
    }
    else {
        return(false);
    }
    
}


long KtreeNode::get_child (char nodeChar) {
    
    int pos = char_to_position(nodeChar);
    
    long child_id = children[pos];
    if (child_id > 0) {
        return(child_id);
    }
    else {
                
        throw("error, no child found from node." + stacktrace());
    }
}

void KtreeNode::add_child (char nodeChar, long id) {

    //cerr << "adding child: " << nodeChar << " : " << id << endl;
    
    int pos = char_to_position(nodeChar);
    
    children[pos] = id;
    
    //cerr << "added child.\n";
}


bool KtreeNode::has_children() {
    
    
    for (int i = 0; i <= 3; i++) {
        if (children[i] > 0) {
            return(true);
        }
    }

    return(false);
}

vector<char> KtreeNode::get_children() {
    
    vector<char> children_chars;

    for (int i = 0; i <= 3; i++) {
        if (children[i] > 0) {
            char nodeChar = position_to_char(i);
            children_chars.push_back(nodeChar);
        }
    }
    
    return(children_chars);
}

    


char KtreeNode::get_char() {
    return(nodeChar);
}


long KtreeNode::get_count() {
    return(this->count);
}

void KtreeNode::set_count(long count) {
    this->count = count;
}


string KtreeNode::toString() {


    stringstream s;

    char c = get_char();
    
    s << "Node[" << c << "]";

    for (int i = 0; i <= 3; i++) {
        long id = children[i];
        if (id > 0) {
        
            char child_c = position_to_char(i);
            s << " :: " << child_c << "," << id;
        }
    }
    s << endl;
    
    return(s.str());
}


int KtreeNode::char_to_position (char node_char) {

    int pos;

    switch (node_char) {

    case 'G':
    case 'g':
        pos = 0;
        break;

    case 'A':
    case 'a':
        pos = 1;
        break;

    case 'T':
    case 't':
        pos = 2;
        break;
        
    case 'C':
    case 'c':
        pos = 3;
        break;

    default:
        stringstream errmsg;
        errmsg << "Error, cannot decipher character: [" << node_char << "] to position. " << stacktrace();
        throw(errmsg.str());
    }


    return(pos);
}
        
char KtreeNode::position_to_char(int pos) {

    char c;

    switch(pos) {

    case 0:
        c = 'G';
        break;
       
    case 1:
        c = 'A';
        break;

    case 2:
        c = 'T';
        break;

    case 3:
        c = 'C';
        break;
        
    default:
        throw("cannot decipher position to a character" + stacktrace());
        
    }

    return(c);
}

