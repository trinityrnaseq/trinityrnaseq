#include "Ktree.hpp"
#include <iostream>
#include <sstream>


Ktree::Ktree() {

    ktree_node_counter = 0;
    
    KtreeNode root_node ('_', 0);
    ktree_node_list.push_back(root_node);

}


KtreeNode& Ktree::get_root_node() {
    
    return(ktree_node_list[0]);
}

void Ktree::add_kmer (string kmer) {
    
    // cerr << "Adding kmer: " << kmer << endl;
    
    long pos = 0; // start with root node.


    for (unsigned int i = 0; i < kmer.length(); i++) {

        KtreeNode& node = ktree_node_list[pos];
                
        char c = kmer[i];
        // cerr << "k-char: " << c << " at pos: " << i << endl;
        // cerr << " corresponds to parent node: " << node.toString();
        
        if (node.has_child(c)) {
            
            pos = node.get_child(c);
            //cerr << "\thas child for: " << c << " at pos: " << pos << ".\n";
            
        }
        else {
            //cerr << "no child, adding it." << endl;
            KtreeNode newNode (c, 0);
            long orig_pos = pos;
            pos = add_node_to_list(newNode); // invalidates earlier node reference, so re-retrieve it below

            KtreeNode& node = ktree_node_list[orig_pos];
            node.add_child(c, pos);

            // cerr << "\tnode now described as: " << node.toString();
            // cerr << "\tOR: " << ktree_node_list[orig_pos].toString();
            
        }
        
        //cerr << "KmerTree:\n" << toString();
        
        //cerr << "reporting kmer counts:" << endl;
        //report_kmer_counts();
    }
    
    // node is now leaf node

    // update the leaf count.
    KtreeNode& node = ktree_node_list[pos];
    node.set_count( node.get_count() + 1);
    ktree_node_list[pos] = node;
    
    // cerr << "-kmer added.\n\n";

}

long Ktree::add_node_to_list(KtreeNode node) {

    // cerr << "Adding node to list." << endl;

    ktree_node_counter++;
    ktree_node_list.push_back(node);
    
    if (ktree_node_list.size() - 1 != ktree_node_counter) {
        throw ("Error, node list and last index are out of synch");
    }

    // cerr << "done.\n";
    return(ktree_node_counter);
    
}

void Ktree::report_kmer_counts() {
    
    KtreeNode& root = get_root_node();
    
    recurse_through_kmer_counts("", root);
}

void Ktree::recurse_through_kmer_counts(string prefix, KtreeNode& node) {

    if (node.has_children()) {
        
        vector<char> children_chars = node.get_children();
        
        for (vector<char>::iterator it = children_chars.begin(); it != children_chars.end(); it++) {
            
            char c = *it;
            long child_id = node.get_child(c);
            KtreeNode& child = ktree_node_list[child_id];
            
            string prefix_extend = prefix + c;
            recurse_through_kmer_counts(prefix_extend, child);
        }
    }

    else {

        cout << prefix << "\t" << node.get_count() << endl;
    }

}

            
            
string Ktree::toString() {

    stringstream s;

    for (vector<KtreeNode>::iterator it = ktree_node_list.begin(); it != ktree_node_list.end(); it++) {
        
        KtreeNode& k = *it;
        
        s << k.toString();
    }

    return(s.str());
}
