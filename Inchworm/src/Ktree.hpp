#ifndef __Ktree__

#define __Ktree__

#include <string>
#include <vector>
#include "KtreeNode.hpp"


using namespace std;

class Ktree {

public:

    Ktree();
    
    
    void add_kmer (string kmer);

    void report_kmer_counts();
    
    string toString();


private:

    
    vector<KtreeNode> ktree_node_list;

    void recurse_through_kmer_counts(string prefix, KtreeNode& node);

    KtreeNode& get_root_node();


    unsigned long ktree_node_counter;
    long add_node_to_list(KtreeNode node);
    
};

#endif

