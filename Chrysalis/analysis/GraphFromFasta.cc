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

static float MIN_KMER_ENTROPY = 1.3;
static float MIN_WELD_ENTROPY = MIN_KMER_ENTROPY;  // min entropy for each half of a welding mer (kk)
static bool DEBUG = false;
static float MAX_RATIO_INTERNALLY_REPETITIVE = 0.85;

static bool REPORT_WELDS = false;
static int MAX_CLUSTER_SIZE = 100;
static int MIN_CONTIG_LENGTH = 24;

static bool __DEBUG_WELD_ALL = false;



// print nucleotide sequence 80 chars per line
void PrintSeq(const DNAVector & d) {
    int i;
    for (i=0; i<d.isize(); i++) {
        cout << d[i];
        if ((i+1) % 80 == 0)
            cout << endl;
    }
    cout << endl;
}


// wrapper around an integer vector
class Pool
{
public:
    
    

    Pool() {
        m_id = -1;
    }
    
    Pool(int id) {
        m_id = id;
    }

    Pool(const Pool& p) {
        m_id = p.m_id;
        m_index = p.m_index;
    }
    
    void push_back(int i) { add(i); }
    void add(int i) {m_index.push_back(i);}
    
    void add (Pool& p) {
        for (int i = 0; i < p.size(); i++) {
            if (! contains(p[i])) {
                add(p[i]);
            }
        }
    }
    
    void clear () {
        svec<int> newvec;
        m_index = newvec;
    }
    
    void exclude (Pool& p) {
        
        for (int i = 0; i < p.size(); i++) {
            if (p.contains(i)) {
                exclude(i);
            }
        }
    }
    
    void exclude (int i) {
        //FIXME:  there must be a more efficient way to do this.
        svec<int> new_mvec;
        for (size_t j = 0; j < m_index.size(); j++) {
            int val = m_index[j];
            if (val != i) {
                new_mvec.push_back(val);
            }
        }
        m_index = new_mvec;
    }
            

    int size() const {return m_index.isize();}
    int isize() const { return(size()); }
    
    int get(int i) const {return m_index[i];}

    int get_id() {
        return(m_id);
    }

    int operator [] (const int& i) {
        return(get(i));
    }

    bool contains(int id) {
        for (int i = 0; i < size(); i++) {
            if (get(i) == id) {
                return(true);
            }
        }
        return(false);
    }
    
    
    void sortvec() {
        sort(m_index.begin(), m_index.end());
    }
    

private:
    svec<int> m_index;
    int m_id;

};

// some prototypes
void add_iworm_link(map<int,Pool>& weld_reinforced_iworm_clusters, int iworm_index, int iworm_pool_addition);
void add_reciprocal_iworm_link(map<int,Pool>& weld_reinforced_iworm_clusters, int iworm_A, int iworm_B);




bool is_simple_repeat(const DNAVector& kmer) {
    
    // in the future, should just do SW alignment
    // for now, compare all half-kmers
    

    string best_left_kmer;
    string best_right_kmer;
    int best_left_pos = 0;
    int best_right_pos = 0;
    float max_repetitive_ratio = 0;

    stringstream left_kmer;
    stringstream right_kmer;

    int mid_kmer_length = kmer.isize()/2;
    
    for (int i = 0; i < mid_kmer_length; i++) {

        for (int j = i+1; j <= mid_kmer_length; j++) {
        
            int ref_kmer_pos = i;
            int internal_kmer_pos = j;
                    

            // compare half-kmers
            
            int bases_compared = 0;
            int bases_common = 0;
     
            left_kmer.clear();
            right_kmer.clear();
            

            while (internal_kmer_pos <= j + mid_kmer_length - 1) {
                bases_compared++;
                if (kmer[ref_kmer_pos] == kmer[internal_kmer_pos]) {
                    bases_common++;
                }
                
                //if (DEBUG) {
                    left_kmer << kmer[ref_kmer_pos];
                    right_kmer << kmer[internal_kmer_pos];
                    //}

                
                ref_kmer_pos++;
                internal_kmer_pos++;
            }
            
            float ratio_same = (float)bases_common/bases_compared;
            if (ratio_same > max_repetitive_ratio) {
                max_repetitive_ratio = ratio_same;
                best_left_kmer = left_kmer.str();
                best_right_kmer = right_kmer.str();
                best_left_pos = i;
                best_right_pos = j;
            }


            if ((! DEBUG) && ratio_same >= MAX_RATIO_INTERNALLY_REPETITIVE) {
                
                // quickest way of exiting this function.  No reason to capture the most repetitive kmer since we're not showing it.
                return(true);
            }
            
        }
           
    }

        
    if (DEBUG) {
        #pragma omp critical
        cout << "# " << kmer.AsString() << " most repetitive kmers:: (" << best_left_pos << "," << best_right_pos << ") " 
             << best_left_kmer << ", " << best_right_kmer << " ratioID: "  << max_repetitive_ratio << endl;
    }

    if (max_repetitive_ratio > MAX_RATIO_INTERNALLY_REPETITIVE) {
        
        if (DEBUG) {
            #pragma omp critical
            cout << "# " << kmer.AsString() << " is internally repetitive: (" << best_left_pos << "," << best_right_pos << ") " 
                 << best_left_kmer << ", " << best_right_kmer << " ratioID: "  << max_repetitive_ratio << endl;
        }
        
        return(true);
    }
    else {
        return(false); // not internally repetitive
    }
}




// used to determine complexity of a kmer
float compute_entropy(const DNAVector & kmer) {
    
    map<char,int> char_map;
    
    for (unsigned int i = 0; i < (unsigned int)kmer.isize(); i++) {
        
        char c = kmer[i];
        char_map[c]++;
    }
    
    float entropy = 0;
    
    char nucs[] = { 'G', 'A', 'T', 'C' };
    
    for (unsigned int i = 0; i < 4; i++) {
        
        char nuc = nucs[i];
        
        int count = char_map[nuc];
        
        float prob = (float)count / kmer.isize();
        
        if (prob > 0) {
            float val = prob * log(1/prob)/log((float)2);
            entropy += val;
        }
    }
    
    return(entropy);
}

bool IsSimple(const DNAVector & d) 
{
   
    // compute entropy
    float e = compute_entropy(d);

    //cout << d.AsString() << " has entropy: " << e << endl;
    if (e < MIN_KMER_ENTROPY) {
        
        //if (DEBUG) 
        //    cerr << "-Sequence declared low complexity: " << d.AsString() << " with entropy " << e << endl;
        
        return(true);
    }
    else {
        return false;
    }
    

    /*    old check
    
    int i;
    int k = 0;
    for (i=2; i<d.isize(); i++) {
        if (d[i] == d[i-2])
            k++;
    }
    double ratio = (double)k/(double)d.isize();
    if (ratio > 0.6) {

        if (DEBUG) 
            cerr << "-Sequence declared low complexity: " << d.AsString() << " with entropy " << e << endl;
    
        return true;
    }

    */
   
    
}


bool SimpleHalves(const DNAVector & d) {
    int len = d.isize();
    int mid_pos = int(len/2);

    DNAVector left;
    left.resize(mid_pos);
    for (int i = 0; i < mid_pos; i++) {
        left[i] = d[i];
    }
    
    DNAVector right;
    right.resize(len - mid_pos);
    int counter = 0;
    for (int i = mid_pos; i < len; i++) {
        right[counter] = d[i];
        
        counter++;
    }

    if (DEBUG) {
        #pragma omp critical
        {
            cout << "## Left half: " << left.AsString() << ", en: " << compute_entropy(left);
            cout << "\tRight half: " << right.AsString() << ", en: " << compute_entropy(right) << endl;
        }
    }
    
    return (compute_entropy(left) < MIN_WELD_ENTROPY 
            || 
            compute_entropy(right) < MIN_WELD_ENTROPY
            ||
            is_simple_repeat(left)
            ||
            is_simple_repeat(right)

            );
}
        


 

// shadows are contigs that are nearly identical but arise from sequencing error-containing kmers
bool IsShadow(const DNAVector & a, const DNAVector & b, int startA, int startB, int k) 
{
    //return false;
    
    
    int i;
    int n = 0;
    int nn = 0;
    int last = -1;
    int len = 0;
    for (i=startA; i<a.isize(); i++) {
        int x = i-startA + startB;
        
        if (x >= b.isize())
            break;
        len++;
        if (a[i] != b[x]) {     
            //cout << "Mismatch @ " << i << " and " << x << endl;
            if (last >= 0) {
                int dist = i-last;
                if (x > 3 && i > 3 && a[i-1] != b[x-1] && a[i-2] != b[x-2])
                    break;
                if (dist == k+1) {
                    n++;
                } else {
                    nn++;
                }
            }
            last = i;
        } else {
            //cout << i << " and " << x << endl;
        }
    }
    
    int expect = (int)(0.9*(double(len/(k+1)-1)));
    //cout << "Len: " << len << " Expect: " << expect << " Observed: " << n << " diss: " << nn << endl;
    if (n >= expect && n > 4 && nn < n/5) {
        return true;
    }
    return false;
    
}


// e.g., given the string:
// >a1;74093 K: 25 length: 8833
// extract "74093"
double Coverage(const string &s) 
{
    const char *nptr = std::strchr(s.c_str(), ';');
    if (nptr == NULL) // ';' not found
        return 1.0;
    
    char *endptr;
    double ret = std::strtod(++nptr, &endptr);
    if (endptr == nptr || ret < 1.0) // conversion not possible, or value < 1.0
        return 1.0;
    
    return ret;
}

bool IsGoodCoverage(double a, double b, double min_iso_ratio) 
{
    
    if (a > b) {
        double tmp = b;
        b = a;
        a = tmp;
    }
    
    if (a/b > min_iso_ratio) {
        return(true);
    }
    else {
        return(false);
    }
        
    /*  original way
    
    //return true;
    
    double dev_a = sqrt(a);
    double dev_b = sqrt(b);
    
    double mean = (a+b)/2;
    
    if (dev_b > dev_a) {
        double dev_tmp = dev_a;
        dev_a = dev_b;
        dev_b = dev_tmp;
        
        double tmp = a;
        a = b;
        b = tmp;
    }
    
    //cout << "Coverage check: " << a << " mean: " << mean << " dev: " << dev_a << endl;
    
    if ((a - mean) < 10*dev_a)
        return true;
    
    double ratio = a/b;
    if (ratio < 1.)
        ratio = 1./ratio;
    
    if (ratio < 100.)
        return true;
    else
        return false;


    */

}



class Welder
{
public:
    Welder(int k, int kk) {
        m_k = k;
        m_kk = kk;
        m_pTab = NULL;
    }
    
    void SetTable(NonRedKmerTable * p) {
        m_pTab = p;
    }
    

    // constructs the required weldable kmer from two contigs a,b and positions one and two
    void WeldableKmer(DNAVector & out, 
                      const DNAVector & a, int one, 
                      const DNAVector & b, int two) 
    {
        out.resize(m_kk);
        int flank = (m_kk - m_k)/2;
        
        int startA = one-flank;
        int stopA = one+m_k;
        int startB = two+m_k;
        int stopB = startB + flank;
        
        if (DEBUG) 
            cerr << "weldableKmer(" << startA << "," << stopA << "); (" << startB << "," << stopB << ")" << endl;
        

        if (startA < 0 || stopB >= b.isize()) {
            out.resize(0);
            if (DEBUG) 
                cerr << "range out of bounds" << endl;
            return;
        }
        
        int i;
        int j = 0;
        for (i=startA; i<stopA; i++) {
            out[j] = a[i];
            j++;
        }
        for (i=startB; i<stopB; i++) {
            out[j] = b[i];
            j++;
        }
        
        if (DEBUG)
            cerr << "\tweld candidate: " << out.AsString() << endl;

    }
    
    
    bool Weldable(const DNAVector & a, int one, const DNAVector & b, int two, int thresh, string& welding_kmer, int& weldable_kmer_read_count) 
    {
        int i;
        DNAVector d; // stores the required weldabler kmer of length kk
        WeldableKmer(d, a, one, b, two); // constructs teh weldable kmer, stores in (d)
        
        if (d.isize() == 0)
            return false; // can't make a weldmer due to out of range.
        

        welding_kmer = d.AsString();
        

        if (IsSimple(d) || SimpleHalves(d)) 
            return false;

        
        int count = m_pTab->GetCount(d, 0); // see if the weldable kmer exists among the reads and get read count.
        weldable_kmer_read_count = count;        
                
        if (thresh == 0)
            return true; // just need a kmer match, dont need read support. (ie. min_glue=0)



        if (DEBUG) 
            cerr << "got weldable candidate: " << d.AsString() << " with count: " << weldable_kmer_read_count << endl;
        

        if (count >= thresh) {

            /* 
               if (SimpleHalves(d)) {
               cerr << "Error, halves of weldmer " << d.AsString() << " were found to fail the simple test....  FIXME" << endl;
               exit(3);
               } 
            */
            
        
            

            return true;
            
        }
        else
            return false;
    }


    int weldmer_count_in_reads(const DNAVector weldmer) 
    {
        // weldmer created outside function and provided as input.
                
        int count = m_pTab->GetCount(weldmer, 0); // see if the weldable kmer exists among the reads and get read count.
        
        return(count);
    }

    
private:
    NonRedKmerTable * m_pTab;
    
    int m_k;
    int m_kk;
};

void Add(vecDNAVector & all, DNAVector & add, int & counter) 
{
    #pragma omp critical
    {
        all.push_back(add);
        counter++;
    }
    
}


void describe_poolings (svec<Pool>& pool_vec) {
    
    cerr << endl << "Pools described as follows:" << endl;
    
    for (unsigned int i=0; i<pool_vec.size(); i++) {
        int oldpool_id = i;
        
        stringstream old_pool_info;
        
        Pool& oldpool = pool_vec[oldpool_id];
        
        old_pool_info << "Pool [" << i << "] info: "; 
        for (int j=0; j<oldpool.size(); j++) {
            int iworm_id = oldpool.get(j);
            
            old_pool_info << iworm_id << " ";
        }
        old_pool_info << endl;
        
        cerr << old_pool_info.str();
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
        
        old_pool_info << endl;
        
        cerr << old_pool_info.str();
    }
    
}


bool sort_pool_sizes_descendingly(const Pool& a, const Pool& b) {
    return (a.size() > b.size());
    
}

bool sort_pool_sizes_ascendingly(const Pool& a, const Pool& b) {
    return (a.size() < b.size());
}


svec<Pool> bubble_up_cluster_growth(map<int,Pool>& pool, map<int,bool>& ) {
    
    vector<Pool> pool_vec;
    for (map<int,Pool>::iterator it = pool.begin(); it != pool.end(); it++) {
        pool_vec.push_back(it->second);
    }
    
    sort (pool_vec.begin(), pool_vec.end(), sort_pool_sizes_ascendingly);
    
    //cerr << "After sorting: " << endl;
    //describe_poolings(pool_vec);
    

    //TODO: nothing ignored ends up as a containment.

    map<int,Pool> pool_idx_to_containment;
    map<int,int>  pool_idx_to_vec_idx;
    // init 
    for (size_t i = 0; i < pool_vec.size(); i++) {
        int id = pool_vec[i].get_id();
        Pool tmp(id);
        pool_idx_to_containment[id] = tmp;
        pool_idx_to_vec_idx[id] = i;
    }
    
    
    cerr << "now bubbling: " << endl << endl;
    bool bubbling = true;
    
    int bubble_round = 0;
    while (bubbling) {
        bubble_round++;
        bubbling = false;
        
        if (DEBUG) {
            describe_bubblings(pool_vec, pool_idx_to_containment, bubble_round);
            cerr << "bubbling up graph, round: " << bubble_round << endl;
        }
        
        for (size_t i = 0; i < pool_vec.size(); i++) {
            
            
            
            Pool& p = pool_vec[i];
            int id = p.get_id();
            
            // cerr << "Bubble_round: " << bubble_round << ", processing pool(" << i << "), with id: " << id << " and size: " << p.size() << endl;
            
            if (pool_idx_to_containment[id].size()) {
                // remove any entries in the pool that are already stored in the containment list.
                p.exclude(pool_idx_to_containment[id]);
            }
            
            if (p.size() > 0) {
                
                // cerr << "Processing pool: " << p.get_id() << " with size: " << p.size() << endl;
                
                // bubble upward
                // get other id:
                bool local_bubbled = false;
                int other_id = -1;
                
                for (int j = 0; j < p.size(); j++) {
                    

                    other_id = p[j];
                    if (other_id == id) { continue; }
                    
                    if (pool_idx_to_containment[other_id].size() + pool_idx_to_containment[id].size() + 2 <= MAX_CLUSTER_SIZE) {  // + 2 since neither self is included in its containment list

                        // move this id to the containment list of the 'other' 
                        pool_idx_to_containment[other_id].add(id);
                        
                        // add this id's containment list to the 'other'
                        pool_idx_to_containment[other_id].add(pool_idx_to_containment[id]);
                        
                        // remove 'other_id' from its own containment list.  (it's self-evident)
                        pool_idx_to_containment[other_id].exclude(other_id);
                        pool_vec[ pool_idx_to_vec_idx[ other_id ] ].exclude( id ); 
                    
                        local_bubbled = true;
                        break;
                    }
                }
                if (local_bubbled) {
                    // update links previously to (id) over to (other_id)

                    for (int j = 0; j < p.size(); j++) {
                        if (p[j] != other_id) {
                            // p[j] should contain id, since id linked to p[j] and all should be reciprocal
                            if (! pool_vec[ pool_idx_to_vec_idx[ p[j] ] ].contains( id )) {
                                

                                describe_bubblings(pool_vec, pool_idx_to_containment, -2);
                                cerr << "Error, " << p[j] << " doesn't contain: " << id << " but vice-versa was true." << endl;
                                exit(4);
                            }
                            
                            pool_vec[ pool_idx_to_vec_idx[ p[j] ] ].exclude( id );
                            // replace it with a link to the other_id, if not already linked.
                            if (! pool_vec[ pool_idx_to_vec_idx[ p[j] ] ].contains( other_id )) {
                                pool_vec[ pool_idx_to_vec_idx[ p[j] ] ].add( other_id );
                            }
                            // links must be reciprocal
                            if (! pool_vec[ pool_idx_to_vec_idx[ other_id ] ].contains( p[j] )) {
                                pool_vec[ pool_idx_to_vec_idx[ other_id ] ].add( p[j] );
                            } 
                            
                        }
                    }
                    
                    // clear out this pool
                    p.clear();
                    
                    // clear out the containment for this pool
                    pool_idx_to_containment[id].clear();
                    
                    
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
                cerr << id << " is a loner, not contained." << endl;
        
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
                
                // cerr << "Mapping: " << " ele: " << ele_id << " to pool " << i << endl;
                
            }
            
        }
    }
    
    int round = 0;
    bool remapped_flag = true;
    while (remapped_flag) {
        


        round++;
        cerr << "Transitive closure, round: " << round << endl;
        remapped_flag = false;
        
        
        if (DEBUG) {
            cerr << "Before remapping: " << endl;
            
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
                    // cerr << "-remapping: iworm(" << iworm_id << ") to pool " << lowest_pool_map_id << endl; 
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
                    //cout << "RELOCATING: " << iworm_id << " from old pool " << oldpool_id << " to new pool " << lowest_pool_map_id << endl;
                }
                //cerr << "-clearing out old pool: " << oldpool_id << endl;
                pool_vec[oldpool_id] = tmp;
                remapped_flag = true;
            }    
        }
        mapped = remapped; // new assignments
        
        if (DEBUG) {
            cerr << "After remapping: " << endl;
            describe_poolings(pool_vec);
        }
        

    }
    
         
    
    //cerr << "...done (" << end_time - start_time << " seconds)" << endl;

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

void add_scaffolds_to_clusters(map<int,Pool>& iworm_clusters, string scaffolding_filename, vecDNAVector& dna, int min_glue_required, float glue_factor) {
    
    ifstream in (scaffolding_filename.c_str());
    
    string line;
    if (in.is_open()) {
        while (! in.eof()) {
            getline(in, line);
        
            if (line.length() == 0) { continue; }
            
            //cerr << line << endl;
            

            /* Format is like so:
            iwormA  idxA  iwormB idxB  count_pair_links

            a12;25 11 a8;4 17 41
            a14;13 13 a3;9 62 40
            a15;15 14 a26;8 25 33
            a13;23 12 a7;4 36 31

            */


            istringstream token (line);
            
            string iworm_acc_A, iworm_index_A_str, iworm_acc_B, iworm_index_B_str, count_str;
            token >> iworm_acc_A >> iworm_index_A_str >> iworm_acc_B >> iworm_index_B_str >> count_str;
            
            int iworm_index_A = atoi(iworm_index_A_str.c_str());
            int iworm_index_B = atoi(iworm_index_B_str.c_str());
            int pair_link_count = atoi(count_str.c_str());
            
            string iworm_A = dna.Name(iworm_index_A);
            string iworm_B = dna.Name(iworm_index_B);

            // verify that our names match up
            if (iworm_A.find(iworm_acc_A) == string::npos) {
                cerr << "Error, cannot locate acc: " << iworm_acc_A << " as substring of " << iworm_A << endl;
                exit(4);
            }
            if (iworm_B.find(iworm_acc_B) == string::npos) {
                cerr << "Error, cannot locate acc: " << iworm_acc_B << " as substring of " << iworm_B << endl;
                exit(4);
            }
            
            
            double coverage_A = Coverage(iworm_A);
            double coverage_B = Coverage(iworm_B);
            
            double higher_coverage_val = (coverage_A > coverage_B) ? coverage_A : coverage_B;
            int minCov = (int) (higher_coverage_val * glue_factor);
            if (minCov < min_glue_required) {
                minCov = min_glue_required;
            }
            
            
            if (pair_link_count >= minCov) {
               
                if (DEBUG || REPORT_WELDS) {
                    cerr << "SCAFFOLD_ACCEPT: " << line << endl;
                    cout << "SCAFFOLD_ACCEPT: " << line << endl;
                }
                
                // reciprocal linkage
                
                add_reciprocal_iworm_link(iworm_clusters, iworm_index_A, iworm_index_B);
                                
            }
            else {
                if (DEBUG) 
                    cerr << "SCAFFOLD_REJECT: " << line << endl;
            }
            
        }
        in.close();
    }
    else {
        cerr << "Error, cannot open file: " << scaffolding_filename;
        exit(3);
    }
    


}


void add_iworm_link(map<int,Pool>& weld_reinforced_iworm_clusters, int iworm_index, int iworm_pool_addition) {
    
    map<int,Pool>::iterator it = weld_reinforced_iworm_clusters.find(iworm_index);

    if (it == weld_reinforced_iworm_clusters.end()) {
        // add it
        Pool p(iworm_index);
        p.add(iworm_pool_addition);
        weld_reinforced_iworm_clusters[iworm_index] = p;
    }
    else {
        Pool& p = it->second;
        if (! p.contains(iworm_pool_addition)) {
            p.add(iworm_pool_addition);
        }
    }
}


void add_reciprocal_iworm_link(map<int,Pool>& weld_reinforced_iworm_clusters, int iworm_A, int iworm_B) {

    add_iworm_link(weld_reinforced_iworm_clusters, iworm_A, iworm_B);
    add_iworm_link(weld_reinforced_iworm_clusters, iworm_B, iworm_A);

}
    

void add_unclustered_iworm_contigs (svec<Pool>& clustered_pools, vecDNAVector& dna) {

    if (DEBUG) {
        cerr << "Adding unclustered iworm contigs." << endl;
    }

    map<int,bool> found;

    for (svec<Pool>::iterator it = clustered_pools.begin(); it != clustered_pools.end(); it++) {
        
        Pool& p = *it;
        
        for (int i = 0; i < p.size(); i++) {
            int id = p[i];
            found[id] = true;
            if (DEBUG)
                cerr << ":found iworm: " << id << " in cluster." << endl;
        }
    }
    

    // add in the missing entries
    for (size_t i = 0; i < dna.size(); i++) {
        if (found.find(i) == found.end()) {
            Pool p;
            p.add(i);
            clustered_pools.push_back(p);
            if (DEBUG)
                cerr << "-iworm: " << i << " added as unclustered." << endl;
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
    // cout << "FATAL ERROR: Could not open file for read: " << s << endl;
    // cout << "Please make sure to enter the correct file name(s). Exiting now." << endl;
    
    return false;
}


int main(int argc,char** argv)
{
    
    
    commandArg<string> aStringCmmd("-i","input fasta");
    commandArg<string> readStringCmmd("-r","read fasta");
    commandArg<bool> strandCmmd("-strand","strand specific", false);
    commandArg<int> kCmmd("-k","k-mer size for pooling", 24);
    commandArg<int> kkCmmd("-kk","k-mer size for welding", 48);
    commandArg<int> threadsCmmd("-t", "number of threads (default: use OMP_NUM_THREADS env var)", -1);
    //commandArg<string> bStringCmmd("-o","output fasta");
    commandArg<string> scaffStringCmmd("-scaffolding", "scaffolded pairs of inchworm contigs", "");
    commandArg<double> glueFactorCmmd("-glue_factor", "fraction of max (iworm pair coverage) for read glue support", 0.05);
    commandArg<int> minGlueCmmd("-min_glue", "absolute min glue support required", 2);
    commandArg<double> minIsoRatioCmmd("-min_iso_ratio", "min ratio of (iworm pair coverage) for join", 0.05);
    commandArg<double> minKmerEntropyCmmd("-min_kmer_entropy", "min entropy value for kmers", MIN_KMER_ENTROPY); 
    commandArg<double> minWeldEntropyCmmd("-min_weld_entropy", "min entropy value for each half of a welding-(kk)mers", MIN_WELD_ENTROPY); 
    commandArg<double> maxRatioInternalRepeatCmmd("-max_ratio_internal_repeat", "maximum ratio identical bases in intra-kmer comparisons", MAX_RATIO_INTERNALLY_REPETITIVE);
    
    commandArg<bool> reportWeldsCmmd("-report_welds", "report the welding kmers", false);
    commandArg<bool> debugCmmd("-debug", "verbosely describes operations", false);
    commandArg<bool> noWeldsCmmd("-no_welds", "no clustering based on weld-reading", false);
    commandArg<int>  maxClusterSizeCmmd("-max_cluster_size", "max size for an inchworm cluster", MAX_CLUSTER_SIZE);
    commandArg<int>  minContigLengthCmmd("-min_contig_length", "min sum cluster contig length", MIN_CONTIG_LENGTH);
    commandArg<bool>  debugWeldAllCmmd("-debug_weld_all", "creates a single cluster of all contigs, for debugging only", false);


    commandLineParser P(argc,argv);
    P.SetDescription("Makes a graph out of a fasta");
    P.registerArg(aStringCmmd);
    P.registerArg(readStringCmmd);
    P.registerArg(strandCmmd);
    P.registerArg(kCmmd);
    P.registerArg(kkCmmd);
    P.registerArg(threadsCmmd);
    P.registerArg(scaffStringCmmd);
    P.registerArg(glueFactorCmmd);
    P.registerArg(minGlueCmmd);
    P.registerArg(minIsoRatioCmmd);
    P.registerArg(maxClusterSizeCmmd);
    
    P.registerArg(minKmerEntropyCmmd);
    P.registerArg(minWeldEntropyCmmd);

    P.registerArg(maxRatioInternalRepeatCmmd);
    P.registerArg(reportWeldsCmmd);
    P.registerArg(noWeldsCmmd);
    P.registerArg(minContigLengthCmmd);
    

    P.registerArg(debugCmmd);
    P.registerArg(debugWeldAllCmmd);

    //P.registerArg(bStringCmmd);
    
    P.parse();
    
    cerr << "-----------------------------------------" << endl
         << "--- Chrysalis: GraphFromFasta -----------" << endl
         << "-- (cluster related inchworm contigs) ---" << endl
         << "-----------------------------------------" << endl << endl;
    


    string iworm_contigs_filename = P.GetStringValueFor(aStringCmmd); //inchworm contigs file
    bool sStrand = P.GetBoolValueFor(strandCmmd); // indicates strand-specific mode
    string readString = P.GetStringValueFor(readStringCmmd); // rna-seq reads file (strand-oriented if in strand-specific mode
    int k = P.GetIntValueFor(kCmmd); // kmer size for pooling, must be 24
    int kk = P.GetIntValueFor(kkCmmd); // kmer size for welding, default 48 = kmer + 1/2 kmer on each side of the kmer
    int num_threads = P.GetIntValueFor(threadsCmmd);
    REPORT_WELDS = P.GetBoolValueFor(reportWeldsCmmd);
    double glue_factor = P.GetDoubleValueFor(glueFactorCmmd);
    int min_glue_required = P.GetIntValueFor(minGlueCmmd);
    double min_iso_ratio = P.GetDoubleValueFor(minIsoRatioCmmd);
    bool bNoWeld = P.GetBoolValueFor(noWeldsCmmd);
    string scaffolding_filename = P.GetStringValueFor(scaffStringCmmd);
    
    MIN_KMER_ENTROPY = P.GetDoubleValueFor(minKmerEntropyCmmd);
    MIN_WELD_ENTROPY = P.GetDoubleValueFor(minWeldEntropyCmmd);
    MAX_RATIO_INTERNALLY_REPETITIVE = P.GetDoubleValueFor(maxRatioInternalRepeatCmmd);
    MAX_CLUSTER_SIZE = P.GetIntValueFor(maxClusterSizeCmmd);
    MIN_CONTIG_LENGTH = P.GetIntValueFor(minContigLengthCmmd);
    
    DEBUG = P.GetBoolValueFor(debugCmmd);
    __DEBUG_WELD_ALL = P.GetBoolValueFor(debugWeldAllCmmd);
    
 
    if (! Exists(readString) ) {
        cerr << "ERROR, cannot open file: " << readString << endl;
        exit(1);
    }
    if (! Exists(iworm_contigs_filename)) {
        cerr << "ERROR, cannot open file: " << iworm_contigs_filename << endl;
        exit(2);
    }
    
    

    if (min_glue_required < 1) {
        cerr << "min_glue set < 1, turning read requirement for welding off: setting min_glue=0, glue_factor=0" << endl;
        min_glue_required = 0;
        glue_factor = 0;
    }
    

    if (num_threads > 0) {
        cerr << "-setting num threads to: " << num_threads << endl;
        omp_set_num_threads(num_threads);
    }

    int real_num_threads;
    #pragma omp parallel
    {
        #pragma omp single
        {
            real_num_threads = omp_get_num_threads();
        }
    }
    cerr << "-running on " << real_num_threads << " threads" << endl;
    
    // read inchworm contigs into memory
    vecDNAVector dna;
    cerr << "GraphFromFasta: Reading file: " << iworm_contigs_filename << endl;
    dna.Read(iworm_contigs_filename, false, false, true, 1000000);
    cerr << "done!" << endl;
    
        
    int i, j; //TODO: remove this, use local vars instead
    

   
    if (k != 24) {
        //cerr << "The only size of k supported is 24! Exiting!" << endl;
        //return -1;
        
        cerr << "Chrysalis::GraphFromFasta requires K=24.... so changing it for just this stage." << endl;
        k = 24; // DEBUG:  WHY 24 required here?  Because it uses sets of 12-mers to define matches, where 2x12mer = 1 24 mer match here.
    }
    
    
    vecDNAVector crossover;
    
    NonRedKmerTable kmers(kk);
    Welder weld(k, kk); // decides if read support exists to weld two inchworm contigs together into the same component.


    int iworm_counter = 0;
    int total_iworm_contigs = dna.size();
    
    
    map<int,Pool> weld_reinforced_iworm_clusters;
    map<int,bool> toasted;
    
    if (scaffolding_filename.length() > 0) {
        // add scaffolding info to clusters
        add_scaffolds_to_clusters(weld_reinforced_iworm_clusters, scaffolding_filename, dna, min_glue_required, glue_factor);
    }

    
                
    // decode inchworm contigs into kmer composition
    KmerAlignCore core;
    core.SetNumTables(2);
    TranslateBasesToNumberExact trans;
    core.SetTranslator(&trans);
    core.AddData(dna);
    core.SortAll();

    //Welder weld(k, kk); // decides if read support exists to weld two inchworm contigs together into the same component.
    

    int schedule_chunksize = dna.size() / real_num_threads / 100;
    if(schedule_chunksize<1) schedule_chunksize=1;
    cerr << "-setting omp for schedule chunksize to " << schedule_chunksize << " for " << dna.size() << " iworm contigs" << endl;
    
    int counter = 0;
    
    double start_time = omp_get_wtime();
    double end_time;

    svec<Pool> clustered_pools;
    
    if (__DEBUG_WELD_ALL) {
    
        // making just one big cluster for testing purposes.
        
        Pool p;
        
        for (size_t i = 0; i < dna.size(); i++) {
            p.add(i);
        }
        
        clustered_pools.push_back(p);
    }

    else {

        if (! bNoWeld) { // make this a function
            cerr << "Phase 1: Collecting candidate weldmers between iworm contig pairs sharing k-1mers" << endl;
            
            iworm_counter = 0;
            
#pragma omp parallel for schedule(dynamic, schedule_chunksize) private(j)
            for (i=0; i<dna.size(); i++) {
                DNAVector & d = dna[i]; // inchworm contig [i]
                
#pragma omp atomic
                iworm_counter++;
                
                if (iworm_counter % 1000 == 0 || iworm_counter == dna.size()-1) {
#pragma omp critical
                    cerr << "\rProcessed: " << iworm_counter/(double)dna.size()*100 << " % of iworm contigs.    ";
                }
                
                for (j=0; j<=d.isize()-k; j++) {
                    DNAVector sub; // a kmer
                    sub.SetToSubOf(d, j, k);
                    
                    if (IsSimple(sub))
                        continue; // ignore kmers that are low complexity
                    
                    
                    // find other inchworm contigs that have kmer matches
                    svec<KmerAlignCoreRecord> matchesFW, matchesRC;   
                    
                    core.GetMatches(matchesFW, sub);
                    if (!sStrand) {
                        sub.ReverseComplement();
                        core.GetMatches(matchesRC, sub);
                    }
                    
                    int x;


                    // create weldmers, and store them for later counting them among the reads.
                    
                    for (x=0; x<matchesFW.isize(); x++) {
                        int c = matchesFW[x].GetContig();
                        if (c == i)
                            continue;  // matches itself
                        
                        DNAVector & dd = dna[c];        
                        int start = matchesFW[x].GetPosition();
                        
                        if (DEBUG) {
                            cerr << "(Phase1: kmer [" << sub.AsString() << "] match between " << dna.Name(i) << " and " << dna.Name(c)
                                 << " at F positions: " << j << " and " << start << endl;
                        } 
                        
                        DNAVector add;
                        
                        weld.WeldableKmer(add, d, j, dd, start);
                        if (add.isize() > 0 && (! IsSimple(add)) && (! SimpleHalves(add) )) {
                            Add(crossover, add, counter);      
                        }
                        weld.WeldableKmer(add, dd, start, d, j);
                        if (add.isize() > 0 && (! IsSimple(add)) && (! SimpleHalves(add) )) {
                            Add(crossover, add, counter);
                        }
                    }
                    for (x=0; x<matchesRC.isize(); x++) {
                        int c = matchesRC[x].GetContig();
                        if (c == i)
                            continue;
                        DNAVector dd = dna[c];
                        dd.ReverseComplement();
                        
                        int start = dd.isize() - matchesRC[x].GetPosition() - k;
                        DNAVector add;
                        
                        weld.WeldableKmer(add, d, j, dd, start); 
                        if (add.isize() > 0 && (! IsSimple(add)) && (! SimpleHalves(add))) {
                            Add(crossover, add, counter);        
                        }
                        
                        
                        weld.WeldableKmer(add, dd, start, d, j);
                        if (add.isize() > 0 && ! IsSimple(add) && (! SimpleHalves(add))) {
                            Add(crossover, add, counter);        
                        }
                    }
                }
            }
            
            end_time = omp_get_wtime();
            cerr << endl << endl << "...done Phase 1. (" << end_time - start_time << " seconds)" << endl;
            
            //crossover.resize(counter);
            
        }
        
        
        
        //----------------------------------------------------
        // -----  Counting the weldmers among the reads ------
        // ---------------------------------------------------
        
        //cerr << "Setting up/sorting structure." << endl;
        cerr << "Captured: " << crossover.size() << " weldmer candidates." << endl;
        
        
        kmers.SetUp(crossover);
        //cerr << "done setting up sorting structure." << endl;
                
        
        cerr << "Setting up reads for streaming..." << endl;
        DNAStringStreamFast seq;
        seq.ReadStream(readString);  //readString is the name of the file containing the reads
        
        cerr << "Identifying reads that support welding of iworm contigs..." << endl;
        kmers.AddData(seq);
        cerr << endl << "Done!" << endl;
        
        weld.SetTable(&kmers);
                
        
        //=================================================================
        // Cluster the inchworm contigs based on read-weld support
        //================================================================
        
        
        cerr << "Phase 2: Reclustering iworm contigs using welds."  << endl;
        
        iworm_counter = 0; // reset
        
#pragma omp parallel for schedule(dynamic, schedule_chunksize) private(j)
        for (i=0; i<dna.size(); i++) {
            
        
#pragma omp atomic
            iworm_counter++;
            
            if (i % 100 == 0) {
#pragma omp critical
                cerr << "\r[" << (iworm_counter/(float)dna.size()*100) << "% done]                ";
            }
            
            
            int cutoff = 0;
            DNAVector & d = dna[i];
            if (d.isize() < k) {
                
                if (DEBUG) {
                    cerr << "ignoring: " << dna.Name(i) << " since less than length: " << k << endl;
                }
                
                continue;
            }
            
            
            
            // iterate through kmers of inchworm contig
            for (j=0; j<=d.isize()-k; j++) {
                
                DNAVector sub;
                sub.SetToSubOf(d, j, k);
                
                svec<KmerAlignCoreRecord> matchesFW, matchesRC;   
                
                core.GetMatches(matchesFW, sub);
                if (!sStrand) {
                    sub.ReverseComplement();
                    core.GetMatches(matchesRC, sub);
                }
                
                if (IsSimple(sub) && matchesFW.isize() + matchesRC.isize() > 1) {
                    if (DEBUG) {
                        cerr << "kmer: " << sub.AsString() << " ignored since either low complex and too many iworm matches"<< endl;
                    }
                    
                    continue;      
                }
                
                int x;
                
                double coverage = Coverage(dna.Name(i));
                
                // Process forward matching kmers
                for (x=0; x<matchesFW.isize(); x++) {
                    int c = matchesFW[x].GetContig();
                    
                    if (c == i) {
                        continue;
                    }
                    
                    
                    if (DEBUG) {
                        cerr << "kmer: " << sub.AsString() << " supports match between " << dna.Name(i) << " and " << dna.Name(c) << endl;
                    }
                    
                    
                    double coverage_other = Coverage(dna.Name(c));
                    
                    
                    double higher_coverage_val = (coverage > coverage_other) ? coverage : coverage_other;
                    int minCov = (int) (higher_coverage_val * glue_factor);
                    if (minCov < min_glue_required) {
                        minCov = min_glue_required;
                    }
                    
                    DNAVector & dd = dna[c];        
                    int start = matchesFW[x].GetPosition();
                    
                    
                    if (DEBUG) {
                        cerr << "Phase2: kmer match [" << sub.AsString() << "] between " << dna.Name(i) << " and " << dna.Name(c)
                             << " at F positions: " << j << " and " << start << endl;
                    } 
                    
                    
                    //cout << "TEST_WELD: " << I << " " << J << " ";
                    string welding_kmer;
                    int welding_kmer_read_count = 0;
                    if ( (! bNoWeld) // welding is on && not weldable, continue
                         && 
                                                  
                         !(weld.Weldable(d, j, dd, start, minCov, welding_kmer, welding_kmer_read_count) 
                           || 
                           weld.Weldable(dd, start, d, j, minCov, welding_kmer, welding_kmer_read_count))) {
                        
                        if (DEBUG) {
                            cerr << "\t no welding kmer avail " << welding_kmer << ", count: " << welding_kmer_read_count 
                                 << ", min count required: " << minCov << endl;
                        }
                        
                        continue;
                    }
                    
                    if (IsShadow(d, dd, j, start, k) && coverage > 2*coverage_other) {
                        cerr << "Toasting shadow: " << dna.Name(c) << endl;
                        
#pragma omp critical
                        toasted[c] = true;
                        
                        continue;
                        
                    } else if (min_glue_required > 0 && !IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
                        //cerr << "Rejecting fw merge between " << dna.Name(i);
                        //cerr << " and " << dna.Name(c) << endl;
                        
                        continue;
                    }
                    // cerr << "Accept (fw)!!" << endl;
                    
#pragma omp critical
                    {
                        
                        add_reciprocal_iworm_link(weld_reinforced_iworm_clusters, i, c);
                        
                        
                    }
                    
                    if (REPORT_WELDS) {
#pragma omp critical
                        cout << "#Welding: " << dna.Name(i) << " to " << dna.Name(c) 
                             << " with kmer:[" << sub.AsString() << "], weldmer:["  << welding_kmer << "] found in " << welding_kmer_read_count << " reads" << endl;
                    }
                    
                    
                }
                
                if (sStrand) {
                    // only doing the forward matches
                    
                    if (DEBUG) {
                        cerr << " only procesing forward strand " << endl;
                    }
                    
                    
                    continue; // superfluous since we didn't capture any rc matches
                }
                
                // Process the RC matches now
                for (x=0; x<matchesRC.isize(); x++) {
                    int c = matchesRC[x].GetContig();
                    
                    
                    if (c == i) {
                        // ignore self matches
                        continue;
                    }
                    
                    
                    if (DEBUG) {
                        cerr << "RCkmer: " << sub.AsString() << " supports match between " << dna.Name(i) << " and " << dna.Name(c) << endl;
                    }
                    
                    
                    
                    double coverage_other = Coverage(dna.Name(c));
                    
                    double higher_coverage_val = (coverage > coverage_other) ? coverage : coverage_other;
                    int minCov = (int) (higher_coverage_val * glue_factor);
                    if (minCov < min_glue_required) {
                        minCov = min_glue_required;
                    }
                    
                    DNAVector dd = dna[c];
                    dd.ReverseComplement();  // revcomp a copy of the iworm[c] sequence
                    
                    int start = dd.isize() - matchesRC[x].GetPosition() - k;  // reverse-complement the match coordinate
                    
                    
                    if (DEBUG) {
                        cerr << "kmer match between " << dna.Name(i) << " and " << dna.Name(c)
                             << " at R positions: " << j << " and " << start << endl;
                    } 
                    
                    
                    string welding_kmer;
                    int welding_kmer_read_count = 0;
                    if (!bNoWeld 
                        && 
                        !(weld.Weldable(d, j, dd, start, minCov, welding_kmer, welding_kmer_read_count) 
                          || 
                          weld.Weldable(dd, start, d, j, minCov, welding_kmer, welding_kmer_read_count))) {
                        
                        
                        if (DEBUG) {
                            cerr << "\t no welding kmer avail " << welding_kmer << ", count: " << welding_kmer_read_count 
                                 << ", min count required: " << minCov << endl;
                        }
                        
                        continue;
                    }
                    
                    if (IsShadow(d, dd, j, start, k) && coverage > 2*coverage_other) {
                        cerr << "Toasting shadow: " << dna.Name(c) << endl;
                        
#pragma omp critical
                        toasted[c] = true;
                        
                        continue;
                        
                    } else if (min_glue_required > 0 && !IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
                        // cerr << "Rejecting rc merge between " << dna.Name(i);
                        // cerr << " and " << dna.Name(c) << endl;
                        continue;
                    }
                    //cerr << "Accept (rc)!!" << endl;
                    
                    
                    
#pragma omp critical
                    {
                        
                        add_reciprocal_iworm_link(weld_reinforced_iworm_clusters, i, c);
                        
                        
                    }
                    
                    if (REPORT_WELDS) {
#pragma omp critical
                        cout << "#Welding: " << dna.Name(i) << " to " << dna.Name(c) 
                             << " with kmer:[" << sub.AsString() << "], weldmer:["  << welding_kmer << "] found in " << welding_kmer_read_count << " reads" << endl;
                    }
                    
                    //cout << "Mapped sequence " << dna.NameClean(c) << " to pool " << mapped[c] << " -" << endl;
                }
                
            }
        }
        
        
        end_time = omp_get_wtime();
        
        cerr << endl; // end of progress monitoring.
        
        
        
        // sl_clustered_pools = sl_cluster_pools(weld_reinforced_iworm_clusters, toasted);
        
        clustered_pools = bubble_up_cluster_growth(weld_reinforced_iworm_clusters, toasted);
        
        if (DEBUG) {
            cerr << "Final pool description: " << endl;
            describe_poolings(clustered_pools);
        }
        
        
        add_unclustered_iworm_contigs(clustered_pools, dna);
                
    }
    
    
    
    //-----------------------------------------------------------------------------------
    // Generate final output
    
    
        
    int component_count = 0;
    
    for (i=0; i<clustered_pools.isize(); i++) {
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

        cout << "COMPONENT " << component_count << "\t" << p.size() << endl;
        for (size_t j = 0; j < p.size(); j++) {
            int z = p.get(j);
            
            pool_info << z << " ";
            cout << ">Component_" << component_count << " " << p.size() << " " << z << " [iworm" << dna.Name(z) << "]" << endl;
            PrintSeq(dna[z]);
        }
        pool_info << endl;
        cout << pool_info.str();
        
        cout << "END" << endl;
    
        component_count++;

    }
    
    return 0;
    
}
  
