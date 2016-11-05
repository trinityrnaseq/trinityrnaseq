#include <string>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <map>
#include <fstream>

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "aligns/KmerAlignCore.h"

#include <math.h>
#include "analysis/KmerTable.h"
#include "analysis/NonRedKmerTable.h"
#include "util/mutil.h"
#include <omp.h>
#include "analysis/DeBruijnGraph.h"
#include "analysis/sequenceUtil.h"
#include "analysis/Pool.h"


static bool DEBUG = false;
static int MAX_CLUSTER_SIZE = 100;
static int MIN_CONTIG_LENGTH = 24;


// print nucleotide sequence 80 chars per line
void PrintSeq(const DNAVector & d) {
    int i;
    for (i=0; i<d.isize(); i++) {
        cout << d[i];
        if ((i+1) % 80 == 0)
            cout << "\n";
    }
    cout << "\n";
}





void Add(vecDNAVector & all, DNAVector & add, int & counter) 
{
    #pragma omp critical
    {
        all.push_back(add);
        counter++;
    }
    
}


void describe_poolings (vector<Pool>& pool_vec) {
    
    cerr << "\n" << "Pools described as follows:" << "\n";
    
    for (unsigned int i=0; i<pool_vec.size(); i++) {
        int oldpool_id = i;
        
        stringstream old_pool_info;
        
        Pool& oldpool = pool_vec[oldpool_id];

        cerr << oldpool.str() << endl;
        
    }
    
}

void describe_bubblings (vector<Pool>& pool_vec, map<int,Pool>& pool_idx_to_containment, int round) {
    
    for (unsigned int i=0; i<pool_vec.size(); i++) {

        stringstream old_pool_info;
        
        Pool& oldpool = pool_vec[i];
        
        old_pool_info << "Bubbling[" << round << "], Pool [" << oldpool.get_id() << "] size: " << oldpool.size() << ", links: "; 
        for (int j=0; j<oldpool.size(); j++) {
            int iworm_id = oldpool.get(j);
            
            old_pool_info << iworm_id << " ";
        }
        if (pool_idx_to_containment.find(oldpool.get_id()) != pool_idx_to_containment.end()) {
            old_pool_info << "; Containments: ";
            Pool& containment_pool = pool_idx_to_containment[oldpool.get_id()];
            for (int j = 0; j < containment_pool.size(); j++) {
                old_pool_info << containment_pool.get(j) << " ";
            }
        }
        
        old_pool_info << "\n";
        
        cerr << old_pool_info.str();
    }
    
}


bool sort_pool_sizes_descendingly(const Pool& a, const Pool& b) {
    return (a.size() > b.size());
    
}

bool sort_pool_sizes_ascendingly(const Pool& a, const Pool& b) {
    return (a.size() < b.size());
}


svec<Pool> bubble_up_cluster_growth(map<int,Pool>& pool) {
    

    /* 
       Algorithm:
       
       A graph is provided as a list of 'node' -> (list of adjacent nodes)
       
       The list is sorted from smallest to largest sets of nodes.
       
       Starting from the smallest entry of {x}  -> { a, b, c}
    
       we start clustering nodes by 'bubbling' up from the smaller clusters to the larger hubs
       by building containment lists and performing node substitutions with containment lists.

       For example, we take 'a' and assign it a containment list [x] and then replace all incoming and outgoing
       edges to x with connections to 'a'.  Now, 'a' represents both 'a' and 'x'.

       We cycle through the lists until there is no more 'bubbling' that is possible.  The 'bubbles' end up representing
       the final cluster memberships.  We put a cap on how big a bubble can grow in order to prevent overclustering.

    */
    
    

    vector<Pool> pool_vec;
    for (map<int,Pool>::iterator it = pool.begin(); it != pool.end(); it++) {
        pool_vec.push_back(it->second);
    }
    
    sort (pool_vec.begin(), pool_vec.end(), sort_pool_sizes_ascendingly);
    
    if (DEBUG) {
        cerr << "After sorting: " << "\n";
        describe_poolings(pool_vec);
    }

    //TODO: nothing ignored ends up as a containment.

    map<int,Pool> pool_idx_to_containment;
    map<int,int>  pool_idx_to_vec_idx;
    // init 
    
    cerr << "bubble_up_cluster_growth(): got " << pool_vec.size() << " pools." << endl;
    
    for (size_t i = 0; i < pool_vec.size(); i++) {
        int pool_id = pool_vec[i].get_id();
        Pool tmp(pool_id);
        pool_idx_to_containment[pool_id] = tmp;
        pool_idx_to_vec_idx[pool_id] = i;
    }
    
    
    cerr << "now bubbling: " << "\n" << "\n";
    bool bubbling = true;
    
    int bubble_round = 0;
    while (bubbling) {
        bubble_round++;
        bubbling = false;
        
        if (DEBUG) {
            describe_bubblings(pool_vec, pool_idx_to_containment, bubble_round);
            cerr << "bubbling up graph, round: " << bubble_round << "\n";
        }
        
        
        // iterate through pools from smaller to larger ones.
        for (size_t i = 0; i < pool_vec.size(); i++) {
                        
            Pool& p = pool_vec[i];
            int id = p.get_id();
            
            // cerr << "Bubble_round: " << bubble_round << ", processing pool(" << i << "), with id: " << id << " and size: " << p.size() << "\n";
            
            if (pool_idx_to_containment[id].size()) {
                // remove any entries in the pool that are already stored in the containment list.
                p.exclude(pool_idx_to_containment[id]);
            }
            
            if (p.size() > 0) {
                
                // cerr << "Processing pool: " << p.get_id() << " with size: " << p.size() << "\n";
                
                // bubble upward
                // get other id:
                bool local_bubbled = false;
                int other_id = -1;
                
                for (int j = 0; j < p.size(); j++) {
                    

                    // find another larger pool we can merge with

                    other_id = p[j];
                    if (other_id == id) {
                        // shouldn't happen since id isn't stored in its own pool.
                        cerr << "Error, other_id: " << other_id << " was stored in its own pool: " << id << endl;
                        continue;
                    } 
                    
                    if (pool_idx_to_containment[other_id].size() + pool_idx_to_containment[id].size() + 2 <= MAX_CLUSTER_SIZE) {  // + 2 since neither self is included in its containment list

                        /////////////////////////////////////////////////////////////////////
                        // begin bubbling of 'id's pool contents to 'other_id's pool contents.
                        /////////////////////////////////////////////////////////////////////

                        if (DEBUG) cerr << "bubbling at id: " << id << ", using " << other_id << " as proxy bubble" << endl;
                        
                        // move this id to the containment list of the 'other' 
                        pool_idx_to_containment[other_id].add(id);
                        
                        // add this id's containment list to the 'other'
                        pool_idx_to_containment[other_id].add(pool_idx_to_containment[id]);
                        
                        // remove 'other_id' from its own containment list.  (it's self-evident)
                        pool_idx_to_containment[other_id].exclude(other_id);

                        // remove 'id' from the other_id's pool
                        pool_vec[ pool_idx_to_vec_idx[ other_id ] ].exclude( id ); 
                        
                        local_bubbled = true;
                        break;
                    }
                }
                if (local_bubbled) {

                    ///////////////////////////////////////////////////////////////////////
                    //  'other_id' is now a proxy for 'id' and all 'id's earlier contents.
                    //   update links previously to (id) over to (other_id)
                    ///////////////////////////////////////////////////////////////////////


                    for (int j = 0; j < p.size(); j++) {
                        
                        int adjacent_id = p[j];
                        Pool& adjacent_pool = pool_vec[ pool_idx_to_vec_idx[ adjacent_id ] ];
                        

                        if (adjacent_id != other_id && adjacent_id != id) {
                            // p[j] should contain id, since id linked to p[j] and all should be reciprocal
                                                        
                            if (! adjacent_pool.contains( id )) {
                                

                                describe_bubblings(pool_vec, pool_idx_to_containment, -2);
                                cerr << "Error, " << p[j] << " doesn't contain: " << id << " but vice-versa was true." << "\n";
                                exit(4);
                            }
                            
                            
                            adjacent_pool.exclude( id );
                            // replace it with a link to the other_id, if not already linked.
                            if (! adjacent_pool.contains( other_id )) {
                                adjacent_pool.add( other_id );
                            }
                            
                            // links must be reciprocal
                            int other_id_pool_vec_idx = pool_idx_to_vec_idx[ other_id ];
                            Pool& other_id_pool = pool_vec[ other_id_pool_vec_idx ];
                            if (! other_id_pool.contains( adjacent_id )) {
                                other_id_pool.add( adjacent_id );
                            } 
                            
                        }
                    } // end of for j iteration through 'id' pool contents.
                    
                    
                    // clear out this pool
                    p.clear();
                    
                    // clear out the containment for this pool since fully covered by other_id's pool now.
                    pool_idx_to_containment[id].clear();


                    // if (DEBUG) describe_bubblings(pool_vec, pool_idx_to_containment, -1);
                                        
                    bubbling = true;
                }
            }
            
        }
        
    }
    
    cerr << "done bubbling.\n";


    // Some entries might be linked but not contained, or exist as singletons.
    map<int,bool> found_contained;
    for (map<int,Pool>::iterator it = pool_idx_to_containment.begin(); it != pool_idx_to_containment.end(); it++) {
        Pool& p = it->second;
        if (p.size() > 0) {
            found_contained[ p.get_id() ] = true; // include self since not in its own containment list.
            for (int i = 0; i < p.size(); i++) {
                int id = p.get(i);
                found_contained[id] = true;
            }
        }
    }

    

    
    svec<Pool> bubbled_up_pools;
    for (size_t i = 0; i < pool_vec.size(); i++) {
                
        Pool& p = pool_vec[i];
        int id = p.get_id();
        
        map<int,Pool>::iterator it = pool_idx_to_containment.find(id);

        if (it != pool_idx_to_containment.end() && (it->second).size() > 0 ) {
            // Entry (id) has a containment list.

            Pool containment_pool = it->second;
            // add (id) and its containment list to the final pool
            containment_pool.add(id);
            bubbled_up_pools.push_back(containment_pool);
        }
        else if (found_contained.find(id) == found_contained.end()) { // not part of any containment list
            
            if (DEBUG) 
                cerr << id << " is a loner, not contained." << "\n";
        
            Pool loner;
            loner.add(id);
            bubbled_up_pools.push_back(loner);
        }
    }

        
    return(bubbled_up_pools);

}

    
// Single-linkage clustering of iworm contigs in pools
svec<Pool> sl_cluster_pools(map<int,Pool>& pool, map<int,bool>& ignore) {

    // Just pull out the ordered pools, initial order doesn't matter... they get trickled upward.
    svec<Pool> pool_vec;
    
    for (map<int,Pool>::iterator it = pool.begin(); it != pool.end(); it++) {
        pool_vec.push_back(it->second);
    }
        
    // init entries to loweset pool index they're found in.
    map<int,int> mapped;
    
    for (size_t i = 0; i < pool_vec.size(); i++) {
        
        Pool p = pool_vec[i];
        
        for (int j = 0; j < p.size(); j++) {
            int ele_id = p.get(j);
            
            if (mapped.find(ele_id) == mapped.end()) {
                
                mapped[ele_id] = i; // first time seeing it, so lowest pool index value
                
                // cerr << "Mapping: " << " ele: " << ele_id << " to pool " << i << "\n";
                
            }
            
        }
    }
    
    int round = 0;
    bool remapped_flag = true;
    while (remapped_flag) {
        


        round++;
        cerr << "Transitive closure, round: " << round << "\n";
        remapped_flag = false;
        
        
        if (DEBUG) {
            cerr << "Before remapping: " << "\n";
            
            describe_poolings(pool_vec);
        }
        

        // do transitive closure on the poolings
        map<int,int> remapped = mapped;
        
        // trickle pool mappings upward
        
        for (unsigned int i=0; i<pool_vec.size(); i++) {
            int oldpool_id = i;
            
            Pool& oldpool = pool_vec[oldpool_id];
            
            int lowest_pool_map_id = oldpool_id;
            bool pool_reassignment_required = false;
            for (int j=0; j<oldpool.size(); j++) {
                int iworm_id = oldpool.get(j);
                int pool_mapping = remapped[iworm_id];
                
                
                if (pool_mapping != lowest_pool_map_id) {
                    pool_reassignment_required = true;
                }
                if (pool_mapping >= 0 && pool_mapping < lowest_pool_map_id) {
                    lowest_pool_map_id = pool_mapping;
                }
            }
            
            // reassign all members to lowest_pool_map_id (if same as pool_id, some entries may still be higher even if not less)
            if (pool_reassignment_required) {
                for (int j=0; j<oldpool.size(); j++) {
                    int iworm_id = oldpool.get(j);
                    remapped[iworm_id] = lowest_pool_map_id;
                    // cerr << "-remapping: iworm(" << iworm_id << ") to pool " << lowest_pool_map_id << "\n"; 
                }
            }
            
            if (lowest_pool_map_id < oldpool_id) {
                // reassign to lower pool
                Pool tmp;
                for (int j=0; j<oldpool.size(); j++) {
                    int iworm_id = oldpool.get(j);
                    
                    if (! pool_vec[lowest_pool_map_id].contains(iworm_id)) {

                        pool_vec[lowest_pool_map_id].add(iworm_id);
                    }
                    //#pragma omp critical
                    //cout << "RELOCATING: " << iworm_id << " from old pool " << oldpool_id << " to new pool " << lowest_pool_map_id << "\n";
                }
                //cerr << "-clearing out old pool: " << oldpool_id << "\n";
                pool_vec[oldpool_id] = tmp;
                remapped_flag = true;
            }    
        }
        mapped = remapped; // new assignments
        
        if (DEBUG) {
            cerr << "After remapping: " << "\n";
            describe_poolings(pool_vec);
        }
        

    }
    
         
    
    //cerr << "...done (" << end_time - start_time << " seconds)" << "\n";

    //------------------------------------------
    // Clustered inchworm contigs now defined.
    // Generate final output
    
    svec<Pool> nr_pools;
    
     for (int i=0; i < (int) pool_vec.size(); i++) {
        Pool & p = pool_vec[i];
        
        if (p.size() == 0) { continue; }
        
        p.sortvec();

        Pool nr_entries; // remove the redundant entries in p
        for (int j=0; j<p.size(); j++) {
            int z = p.get(j);
            
            if (ignore[z]) {
                continue;
            }
            

            if (j > 0 && p.get(j-1) == z) {
                // same entry ended up on the pool vec, already reported it.
                continue;
            }
            
            nr_entries.add(z);
        }
        nr_pools.push_back(nr_entries);
    }
    

    return(nr_pools);
    

}


void add_unclustered_iworm_contigs (svec<Pool>& clustered_pools, vecDNAVector& dna) {

    if (DEBUG) {
        cerr << "Adding unclustered iworm contigs." << "\n";
    }

    map<int,bool> found;

    for (svec<Pool>::iterator it = clustered_pools.begin(); it != clustered_pools.end(); it++) {

        Pool& p = *it;

        for (int i = 0; i < p.size(); i++) {
            int id = p[i];
            found[id] = true;
            if (DEBUG)
                cerr << ":found iworm: " << id << " in cluster." << "\n";
        }
    }


    // add in the missing entries
    for (size_t i = 0; i < dna.size(); i++) {
        if (found.find(i) == found.end()) {
            Pool p;
            p.add(i);
            clustered_pools.push_back(p);
            if (DEBUG)
                cerr << "-iworm: " << i << " added as unclustered." << "\n";
        }
    }

}



bool Exists(const string & s) 
{
    FILE * p = fopen(s.c_str(), "r");
    if (p != NULL) {
        fclose(p);
        return true;
    }
    // cout << "FATAL ERROR: Could not open file for read: " << s << "\n";
    // cout << "Please make sure to enter the correct file name(s). Exiting now." << "\n";
    
    return false;
}


void populate_weld_reinforced_iworm_clusters(string& weld_graph_file, map<int,Pool>& weld_reinforced_iworm_clusters) {

    cerr << "Parsing weld graph file: " << weld_graph_file << endl;
    
    ifstream in (weld_graph_file.c_str());
    if (! in.is_open()) {
        cerr << "Error, cannot open file: " << weld_graph_file << endl;
        exit(3);
    }

    // lines look like this:
    //    210 -> 735
    //    907 -> 242 506
    //    735 -> 515 13 807 779 85 92 210 372 397 567 819    
    
    while (! in.eof()) {
        string line;
        getline(in, line);
        //cerr << "input_line: " << line << "\n";
        
        istringstream token (line);

        string tok;
        token >> tok;

        int node_id = atoi(tok.c_str());
        Pool p(node_id);
        
        token >> tok; // pointer '->' separator string.

        
        while (! token.eof()) {
            token >> tok;
            int adjacent_node = atoi(tok.c_str());
            p.add(adjacent_node);
            //cerr << "Tok: [" << tok << "]" << "\n";
        }
        
        if (node_id >= 0 && p.size() > 0) {
            weld_reinforced_iworm_clusters[ node_id ] = p;
            if (DEBUG) {
                cerr << "Assigning node_id [" << node_id << "] to pool: " << p.str() << endl;
            }
        }
    }
    
    
    
    return;
}



int main(int argc,char** argv)
{
    
    
    commandArg<string> aStringCmmd("-i","input fasta");
    commandArg<string> weldGraphCmmd("-weld_graph", "iworm index weld graph");
    commandArg<bool> debugCmmd("-debug", "verbosely describes operations", false);
    commandArg<int>  minContigLengthCmmd("-min_contig_length", "min sum cluster contig length", MIN_CONTIG_LENGTH);
    
    
    commandLineParser P(argc,argv);
    P.SetDescription("Makes a graph out of a fasta");
    P.registerArg(aStringCmmd);
    P.registerArg(weldGraphCmmd);
    P.registerArg(debugCmmd);
    P.registerArg(minContigLengthCmmd);
    
    P.parse();
    
    cerr << "-------------------------------------------------------" << "\n"
         << "--- Chrysalis: BubbleUpClustering ---------------------" << "\n"
         << "-- (generating the final clusters of iworm contigs) ---" << "\n"
         << "-------------------------------------------------------" << "\n" << "\n";
    
    string iworm_contigs_filename = P.GetStringValueFor(aStringCmmd); //inchworm contigs file
    DEBUG = P.GetBoolValueFor(debugCmmd);
    string weld_graph_file = P.GetStringValueFor(weldGraphCmmd);
    MIN_CONTIG_LENGTH = P.GetIntValueFor(minContigLengthCmmd);
    
    
    if (! Exists(iworm_contigs_filename)) {
        cerr << "ERROR, cannot open iworm contigs file: " << iworm_contigs_filename << "\n";
        exit(2);
    }
    if (! Exists(weld_graph_file)) {
        cerr << "ERROR, cannot open weld graph file: " << weld_graph_file << "\n";
        exit(3);
    }
    
    // read inchworm contigs into memory
    vecDNAVector dna;
    cerr << "BubbleUpClustering: Reading file: " << iworm_contigs_filename << "\n";
    dna.Read(iworm_contigs_filename, false, false, true, 1000000);
    cerr << "done!" << "\n";
    
        
    map<int,Pool> weld_reinforced_iworm_clusters;

    populate_weld_reinforced_iworm_clusters(weld_graph_file, weld_reinforced_iworm_clusters);
    
    svec<Pool> clustered_pools = bubble_up_cluster_growth(weld_reinforced_iworm_clusters);
        
    if (DEBUG) {
        cerr << "Final pool description: " << "\n";
        describe_poolings(clustered_pools);
    }

        
    add_unclustered_iworm_contigs(clustered_pools, dna);
    
        
    //-----------------------------------------------------------------------------------
    // Generate final output
    
    
        
    int component_count = 0;
    
    for (int i=0; i<clustered_pools.isize(); i++) {
        Pool & p = clustered_pools[i];
        
        if (p.size() == 0) { continue; }

        p.sortvec();
    
        int sum_iworm_length = 0;
        for (size_t j = 0; j < p.size(); j++) {
            int z = p.get(j);
            sum_iworm_length += dna[z].isize();
        }

        if (sum_iworm_length < MIN_CONTIG_LENGTH) {
            continue;
        }
        
        stringstream pool_info;                
        pool_info << "#POOL_INFO\t" << component_count << ":" << "\t";

        cout << "COMPONENT " << component_count << "\t" << p.size() << "\n";
        for (size_t j = 0; j < p.size(); j++) {
            int z = p.get(j);
            
            pool_info << z << " ";
            cout << ">Component_" << component_count << " " << p.size() << " " << z << " [iworm" << dna.Name(z) << "]" << "\n";
            PrintSeq(dna[z]);
        }
        pool_info << "\n";
        cout << pool_info.str();
        
        cout << "END" << "\n";
    
        component_count++;

    }
    
    return 0;
    
}
  
