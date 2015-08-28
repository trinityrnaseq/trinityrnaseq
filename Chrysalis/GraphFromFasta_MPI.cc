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

#define HUGE_DATASET
#ifdef MPI_ENABLED
#include<mpi.h>
#endif
#include <unistd.h>
#include <sys/resource.h>
void print_cpu_memory()
{
 struct rusage usage;
 getrusage (RUSAGE_SELF, &usage);
 // printf ("CPU time: %ld.%06ld sec user, %ld.%06ld sec system\n",
 //  usage.ru_utime.tv_sec, usage.ru_utime.tv_usec,
 //    usage.ru_stime.tv_sec, usage.ru_stime.tv_usec);

 cerr<<"Max memory consumed  in GraphFromFasta "<<usage.ru_maxrss<<endl;
 

}



static float MIN_KMER_ENTROPY = 1.3;
static float MIN_WELD_ENTROPY = MIN_KMER_ENTROPY;  // min entropy for each half of a welding mer (kk)
static bool DEBUG = false;
static float MAX_RATIO_INTERNALLY_REPETITIVE = 0.85;

static bool REPORT_WELDS = false;
static int MAX_CLUSTER_SIZE = 100;
static int MIN_CONTIG_LENGTH = 24;


#include<sys/time.h>
static struct timeval start,end;
void timer_start(){

  gettimeofday(&start,NULL);
}


double timer_stop(){
  gettimeofday(&end,NULL);
  double time_taken = ((end.tv_usec-start.tv_usec) + 1000000*(end.tv_sec - start.tv_sec));
 time_taken = time_taken/1000;
  return time_taken;
}


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
        for (int j = 0; j < m_index.size(); j++) {
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


    int mid_kmer_length = kmer.isize()/2;
    
    for (int i = 0; i < mid_kmer_length; i++) {

        for (int j = i+1; j <= mid_kmer_length; j++) {
        
            int ref_kmer_pos = i;
            int internal_kmer_pos = j;
                    

            // compare half-kmers
            
            int bases_compared = 0;
            int bases_common = 0;
     
            stringstream left_kmer;
            stringstream right_kmer;
            

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
            return false;
        
        
        int count = m_pTab->GetCount(d, 0); // see if the weldable kmer exists among the reads and get read count.
        weldable_kmer_read_count = count;
        

        if (DEBUG) 
            cerr << "got weldable candidate: " << d.AsString() << " with count: " << weldable_kmer_read_count << endl;
        

        if (count >= thresh) {

            /* 
               if (SimpleHalves(d)) {
               cerr << "Error, halves of weldmer " << d.AsString() << " were found to fail the simple test....  FIXME" << endl;
               exit(3);
               } 
            */
            
            welding_kmer = d.AsString();

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


svec<Pool> bubble_up_cluster_growth(map<int,Pool>& pool, map<int,bool>& ignore) {
    
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
    for (int i = 0; i < pool_vec.size(); i++) {
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
        
        for (int i = 0; i < pool_vec.size(); i++) {
            
            
            
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
    for (int i = 0; i < pool_vec.size(); i++) {
                
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
               
                if (DEBUG)
                    cerr << "SCAFFOLD_ACCEPT: " << line << endl;
            
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
        
        for (size_t i = 0; i < p.size(); i++) {
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


int main(int argc,char** argv)
{
  int rank,numranks;
#ifdef MPI_ENABLED
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
  rank=0;
  numranks=1;
#endif


 
    
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
    commandArg<bool> noWeldsCmmd("-no_welds", "disable requirement for glue support", false);
    commandArg<int>  maxClusterSizeCmmd("-max_cluster_size", "max size for an inchworm cluster", MAX_CLUSTER_SIZE);
    commandArg<int>  minContigLengthCmmd("-min_contig_length", "min sum cluster contig length", MIN_CONTIG_LENGTH);

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
    
    //P.registerArg(bStringCmmd);
    
    P.parse();
    
    cerr << "-----------------------------------------" << endl
         << "--- Chrysalis: GraphFromFasta -----------" << endl
         << "-- (cluster related inchworm contigs) ---" << endl
         << "-----------------------------------------" << endl << endl;
    


    string aString = P.GetStringValueFor(aStringCmmd); //inchworm contigs file
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

    if (min_glue_required < 1) {
        cerr << "-error, cannot have less than 1 read as glue.  Setting -min_glue to 1" << endl;
        min_glue_required = 1;
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
    cerr << "GraphFromFasta: Reading file: " << aString << endl;
    dna.Read(aString, false, false, true, 1000000);
    cerr << "done!" << endl;
    
        
    size_t i, j; //TODO: remove this, use local vars instead
    
    if (k != 24) {
        cerr << "The only size of k supported is 24! Exiting!" << endl;
        return -1;
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

#ifdef MPI_ENABLED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (! bNoWeld) { // make this a function
        
        
                
        // decode inchworm contigs into kmer composition
        KmerAlignCore core;
        core.SetNumTables(2);
        TranslateBasesToNumberExact trans;
        core.SetTranslator(&trans);
        core.AddData(dna);
        core.SortAll();
            
        cerr << "Phase 1: Collecting candidate weldmers between iworm contig pairs sharing k-1mers" << endl;
        

        Welder weld(k, kk); // decides if read support exists to weld two inchworm contigs together into the same component.
        
#ifdef MPI_ENABLED
        int schedule_chunksize = dna.size() / (numranks) / 100;
#else
	int schedule_chunksize = dna.size() / (real_num_threads) / 100; 
#endif


        if(schedule_chunksize<1) schedule_chunksize=1;
        cerr << "-setting omp for schedule chunksize to " << schedule_chunksize << " for " << dna.size() << " iworm contigs" << endl;

        int counter = 0;
        
        double start_time = omp_get_wtime();
        
        iworm_counter = 0;


	int rank_chunksize;
	
	timer_start();

	int thread_chunksize=schedule_chunksize /(real_num_threads) / 10;
	int thread_count;

#ifdef MPI_ENABLED
	 for (int mpi_count=rank*schedule_chunksize; mpi_count<dna.isize(); mpi_count+=numranks*schedule_chunksize) {
	  if( (dna.isize()-mpi_count)<schedule_chunksize)
	    rank_chunksize=dna.isize()-mpi_count;
	  else
	    rank_chunksize=schedule_chunksize;

#pragma omp parallel for schedule(static) private(i,j,thread_count)
	  for(thread_count=0;thread_count<rank_chunksize;thread_count++){

	    i=mpi_count + thread_count;


	    DNAVector & d = dna[i]; // inchworm contig [i]
            
            #pragma omp atomic
            iworm_counter++;

            if (iworm_counter % 1000 == 0 || iworm_counter == dna.isize()-1) {
            #pragma omp critical
	      cerr << "\rProcessed: " << iworm_counter/(double)dna.isize()*100 << " % of iworm contigs on rank "<<rank   ;
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
            
                for (x=0; x<matchesFW.isize(); x++) {
                    int c = matchesFW[x].GetContig();
                    if (c == i)
                        continue;  // matches itself
                
                    DNAVector & dd = dna[c];        
                    int start = matchesFW[x].GetPosition();
                
                    if (DEBUG) {
                        cerr << "(Phase1: kmer match between " << dna.Name(i) << " and " << dna.Name(c)
                             << " at F positions: " << j << " and " << start << endl;
                    } 

                    DNAVector add;
                
                    weld.WeldableKmer(add, d, j, dd, start);
                    if (add.isize() > 0 && ! IsSimple(add) && (! SimpleHalves(add) )) {
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
                    if (add.isize() > 0 && (! IsSimple(add)) && (! SimpleHalves(add)))
                        Add(crossover, add, counter);        
                
                
                    weld.WeldableKmer(add, dd, start, d, j);
                    if (add.isize() > 0 && ! IsSimple(add) && (! SimpleHalves(add)))
                        Add(crossover, add, counter);        
                }
            }
	  }
	 }
#else
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
            
                for (x=0; x<matchesFW.isize(); x++) {
                    int c = matchesFW[x].GetContig();
                    if (c == i)
                        continue;  // matches itself
                
                    DNAVector & dd = dna[c];        
                    int start = matchesFW[x].GetPosition();
                
                    if (DEBUG) {
                        cerr << "(Phase1: kmer match between " << dna.Name(i) << " and " << dna.Name(c)
                             << " at F positions: " << j << " and " << start << endl;
                    } 

                    DNAVector add;
                
                    weld.WeldableKmer(add, d, j, dd, start);
                    if (add.isize() > 0 && ! IsSimple(add) && (! SimpleHalves(add) )) {
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
                    if (add.isize() > 0 && (! IsSimple(add)) && (! SimpleHalves(add)))
                        Add(crossover, add, counter);        
                
                
                    weld.WeldableKmer(add, dd, start, d, j);
                    if (add.isize() > 0 && ! IsSimple(add) && (! SimpleHalves(add)))
                        Add(crossover, add, counter);        
                }
            }
        }


#endif        
	 
	 
        double end_time = omp_get_wtime();
        cerr << endl << endl << "...done Phase 1. (" << end_time - start_time << " seconds)" << endl;

	cerr<< endl << endl << "done Phase 1 with my timers "<<timer_stop()<<" on rank "<<rank<<" (in milliseconds)"<<endl;
	cerr<<"Crossover size on rank "<<rank<<" is "<<crossover.size()<<std::endl;
        
        //crossover.resize(counter);
        
        
	timer_start();
#ifdef MPI_ENABLED
	


	char* combined_dnas_chars;

#ifdef HUGE_DATASET
	long* dnachars_counts, *crossover_counts;
	dnachars_counts = (long*)malloc(sizeof(long)*numranks);
	long crossover_size = crossover.isize();
	crossover_counts = (long*)malloc(sizeof(long)*numranks);
#else
	int* dnachars_counts, *crossover_counts;
	dnachars_counts = (int*)malloc(sizeof(int)*numranks);
	int crossover_size = crossover.isize();
	crossover_counts = (int*)malloc(sizeof(int)*numranks);
#endif


	  


	

	  



	 
	
	  DNAVector combined_dnavector;
	  for(int count=0;count<crossover.isize();count++)
	      combined_dnavector +=crossover[count];

#ifdef HUGE_DATASET
	  long dnachars_size = combined_dnavector.size();
#else
	  int dnachars_size = combined_dnavector.size();
#endif
	

	  cerr<<std::endl<<"Total dna chars  "<<combined_dnavector.size()<<" on rank "<<rank<<" with a total memory of "<<((double)(combined_dnavector.size()*sizeof(char))/(double)(1024*1024))<<" (MB)"<<std::endl;

	  char* dnas_chars_send;
	  svec<char>dnadata = combined_dnavector.Data();

	  dnas_chars_send = (char*)malloc(sizeof(char)*combined_dnavector.size());
	  for(long count=0;count<combined_dnavector.size();count++)
	    dnas_chars_send[count] = dnadata[count];
	    



	
#ifdef HUGE_DATASET
	  MPI_Allgather(&crossover_size,1,MPI_LONG,crossover_counts,1,MPI_LONG,MPI_COMM_WORLD);
          MPI_Allgather(&dnachars_size, 1,MPI_LONG,dnachars_counts, 1,MPI_LONG,MPI_COMM_WORLD);	  
#else
	  MPI_Allgather(&crossover_size,1,MPI_INT,crossover_counts,1,MPI_INT,MPI_COMM_WORLD);
          MPI_Allgather(&dnachars_size, 1,MPI_INT,dnachars_counts, 1,MPI_INT,MPI_COMM_WORLD);

#endif


	  

	
	  long long combined_crossover_counts=0;

	    for(int count=0;count<numranks;count++)
	      combined_crossover_counts += crossover_counts[count];
	  
	    long long combined_dnachars_counts=0;

	    for(int count=0;count<numranks;count++)
	      {
		if(rank==0)
		  cerr<<"DNA["<<count<<"]="<<dnachars_counts<<endl;
		combined_dnachars_counts += dnachars_counts[count];

	      }


	
	    cerr<<std::endl<<"Combined DNA characters to be exchanged "<<combined_dnachars_counts<<" on rank "<<rank<<" with a total memory of "<<((double)(combined_dnachars_counts*sizeof(char))/(double)(1024*1024))<<" (MB)"<<std::endl;

	    MPI_Barrier(MPI_COMM_WORLD);
	    combined_dnas_chars = (char*)malloc(sizeof(char)*combined_dnachars_counts);
	   


#ifdef HUGE_DATASET
	       long* dnachar_displs;
	       dnachar_displs = (long*)malloc(sizeof(long)*numranks);
#else
	       int *dnachar_displs;
	       dnachar_displs = (int*)malloc(sizeof(int)*numranks);
#endif


	    
	      dnachar_displs[0]=0;
	       for(int count=1;count<numranks;count++)
		 dnachar_displs[count] = dnachar_displs[count-1]+dnachars_counts[count-1];

	       cerr<<"Allocated memory for all DNA characters on rank "<<rank<<std::endl;

	       
#ifdef HUGE_DATASET
	    
	     
	       for(long count=0;count<dnachars_counts[rank];count++)
		 combined_dnas_chars[dnachar_displs[rank]+count]=dnas_chars_send[count];
	       
	       for(int count=0;count<numranks;count++)
		 {
		   void* buffer = &combined_dnas_chars[dnachar_displs[count]];
		   MPI_Bcast(buffer,dnachars_counts[count],MPI_CHAR,count,MPI_COMM_WORLD);
		   
		 }
	       

#else
	       MPI_Allgatherv(dnas_chars_send, dnachars_size, MPI_CHAR,  combined_dnas_chars, dnachars_counts, dnachar_displs,MPI_CHAR,MPI_COMM_WORLD);
#endif

	  free(dnas_chars_send);
	  

	  vecDNAVector combined_crossover;

	      int t=0;
	       for(long long count=0;count<combined_crossover_counts;count++)
		 { 

		   const string s(&combined_dnas_chars[t],kk);

		   //	   cerr<<"string = "<<s<<std::endl;
		   DNAVector dna_sequence;
		   dna_sequence.SetFromBases(s);
		   combined_crossover.push_back(dna_sequence);

		   t = t + kk;

		 }
	      
	       cerr<<"Subsequence on rank "<<rank<<" in beginning of combined DNA characters (8) = ";
	       for(int count=0;count<8;count++)
		 cerr<<combined_dnas_chars[count];

	       cerr<<std::endl;

	    
	        cerr<<"Subsequence on rank "<<rank<<" in end of combined DNA characters (8) = ";
	       for(int count=0;count<8;count++)
		 cerr<<combined_dnas_chars[combined_dnachars_counts-count];

	       cerr<<std::endl;


	       cerr<<"Size of combined crossover = "<<combined_crossover.isize()<<endl;
	       free(combined_dnas_chars);
	       free(crossover_counts);
	       free(dnachar_displs);
	       free(dnachars_counts);


	       
	       
	



#endif

	  if(rank==0)
	    cerr<< endl << endl << " Overhead of MPI steps  "<<timer_stop()<<" (in milliseconds)"<<endl;
 
        //-------------------------------------------------------------------------
        //------------  Now, do iworm clustering using welds ----------------------
        //-------------------------------------------------------------------------
#ifdef MPI_ENABLED        
        kmers.SetUp(combined_crossover);
#else
        kmers.SetUp(crossover);
#endif


        //cerr << "done setting up sorting structure." << endl;
    
        cerr << "Setting up reads for streaming..." << endl;
        DNAStringStreamFast seq;
        seq.ReadStream(readString);  //readString is the name of the file containing the reads
        
        cerr << "Identifying reads that support welding of iworm contigs..."<<endl;
	
        kmers.AddData(seq);
        cerr << endl << "Done!" << endl;
        
        weld.SetTable(&kmers);

#ifdef MPI_ENABLED
        MPI_Barrier(MPI_COMM_WORLD);
#endif

	start_time = omp_get_wtime();
        timer_start();
    
        //=================================================================
        // Cluster the inchworm contigs based on read-weld support
if(rank==0)
        cerr << "Phase 2: Reclustering iworm contigs using welds."  << endl;

#ifdef MPI_ENABLED

	svec<int> ind_leftFW;
	svec<int> ind_rightFW;
	svec<int> ind_leftRC;
	svec<int> ind_rightRC;	
	svec<int> ind_toasted;

	iworm_counter = 0; // reset 

         for (int mpi_count=rank*schedule_chunksize; mpi_count<dna.isize(); mpi_count+=numranks*schedule_chunksize) {
          if( (dna.isize()-mpi_count)<schedule_chunksize)
            rank_chunksize=dna.isize()-mpi_count;
          else
            rank_chunksize=schedule_chunksize;

#pragma omp parallel for schedule(static) private(j, thread_count,i)
          for(thread_count=0;thread_count<rank_chunksize;thread_count++) {

            i=mpi_count + thread_count;

	    #pragma omp atomic
            iworm_counter++;

            if (i % 100 == 0) {
                #pragma omp critical
                cerr << "\r[" << (iworm_counter/(float)dna.isize()*100) << "% done]                ";
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
                        cerr << "Phase2: kmer match between " << dna.Name(i) << " and " << dna.Name(c)
                             << " at F positions: " << j << " and " << start << endl;
                    }


		    //cout << "TEST_WELD: " << I << " " << J << " ";
		     string welding_kmer;
                    int welding_kmer_read_count;
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
                    //    toasted[c] = true;
                          ind_toasted.push_back(c);

                        continue;

                    } else if (!IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
			//cerr << "Rejecting fw merge between " << dna.Name(i);
			//cerr << " and " << dna.Name(c) << endl;
			
			 continue;
                    }
		    // cerr << "Accept (fw)!!" << endl;

		    #pragma omp critical
                    {

                      //  add_reciprocal_iworm_link(weld_reinforced_iworm_clusters, i, c);
			ind_leftFW.push_back(i);
			ind_rightFW.push_back(c);
                    }

                    if (REPORT_WELDS) {
                    #pragma omp critical
                        cout << "#Welding: " << dna.Name(i) << " to " << dna.Name(c)
                             << " with " << welding_kmer << " found in " << welding_kmer_read_count << " reads" << endl;
                    }


                }

                if (sStrand) {
			// only doing the forward matches

		    if (DEBUG) {
                        cerr << " only procesing forward strand " << endl;
                    }

                    continue;   // superfluous since we didn't capture any rc matches
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
                    dd.ReverseComplement(); // revcomp a copy of the iworm[c] sequenc

		    int start = dd.isize() - matchesRC[x].GetPosition() - k; // reverse-complement the match coordinate

		    if (DEBUG) {
                        cerr << "kmer match between " << dna.Name(i) << " and " << dna.Name(c)
                             << " at R positions: " << j << " and " << start << endl;
                    }


                    string welding_kmer;
                    int welding_kmer_read_count;
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
                        // toasted[c] = true;
			ind_toasted.push_back(c);

                        continue;

                    } else if (!IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
			// cerr << "Rejecting rc merge between " << dna.Name(i);
			// cerr << " and " << dna.Name(c) << endl;
			continue;
                    } else if (!IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
			// cerr << "Rejecting rc merge between " << dna.Name(i);
			// cerr << " and " << dna.Name(c) << endl;
			continue;
                    }
			//cerr << "Accept (rc)!!" << endl;


		   #pragma omp critical
                    {

                      //  add_reciprocal_iworm_link(weld_reinforced_iworm_clusters, i, c);
			ind_leftRC.push_back(i);
			ind_rightRC.push_back(c);
                    }

                    if (REPORT_WELDS) {
                    #pragma omp critical
                        cout << "#Welding: " << dna.Name(i) << " to " << dna.Name(c)
                             << " with " << welding_kmer << " found in " << welding_kmer_read_count << " reads" << endl;
                    }

		   //cout << "Mapped sequence " << dna.NameClean(c) << " to pool " << mapped[c] << " -" << endl;
	       }

            }

        }
	}


	end_time = omp_get_wtime();
	cerr << endl << endl << "...done Phase 2 on rank "<< rank << " (" << end_time - start_time << " seconds)" << endl;

	cerr<<"\nPhase 2 took "<<timer_stop()<<" on rank "<<rank<<" (in milliseconds)"<<endl;
        cerr << endl; // end of progress monitoring.


        MPI_Barrier(MPI_COMM_WORLD); 
  	timer_start();

	int* all_indexFW_sizes;
	int* all_indexRC_sizes;
	int* all_IndexFW_left;
	int* all_IndexFW_right;
	int* all_IndexRC_left;
	int* all_IndexRC_right;
	int* leftFW;
	int* rightFW;
	int* leftRC;
	int* rightRC;
	int numFW, numRC, numToasted;
	int* Ind_Toasted;
	int* all_Toasted_sizes;
	int* all_Index_Toasted;

	numToasted = (int)ind_toasted.size();
	Ind_Toasted = (int*)malloc(sizeof(int)*numToasted);
	for(int count=0; count<numToasted; count++)
		Ind_Toasted[count] = ind_toasted[count];	

	numFW = (int)ind_leftFW.size();
	if(!sStrand)	
	   numRC = (int)ind_leftRC.size();

	leftFW = (int*)malloc(sizeof(int)*numFW);
	rightFW = (int*)malloc(sizeof(int)*numFW);
	for(int count=0; count<numFW; count++) {
		leftFW[count]  = ind_leftFW[count];
		rightFW[count] = ind_rightFW[count]; 
	}

	if (!sStrand) {
	  leftRC = (int*)malloc(sizeof(int)*numRC);
          rightRC = (int*)malloc(sizeof(int)*numRC);
          for(int count=0; count<numRC; count++) {
                leftRC[count]  = ind_leftRC[count];
                rightRC[count] = ind_rightRC[count];
       	  } 
	}


	all_indexFW_sizes = (int*)malloc(sizeof(int)*numranks); 
	if (!sStrand) 
	  all_indexRC_sizes = (int*)malloc(sizeof(int)*numranks);

	all_Toasted_sizes =  (int*)malloc(sizeof(int)*numranks);


	MPI_Allgather(&numFW,1,MPI_INT,all_indexFW_sizes,1,MPI_INT,MPI_COMM_WORLD);
	if (!sStrand)
	   MPI_Allgather(&numRC,1,MPI_INT,all_indexRC_sizes,1,MPI_INT,MPI_COMM_WORLD);

	MPI_Allgather(&numToasted,1,MPI_INT,all_Toasted_sizes,1,MPI_INT,MPI_COMM_WORLD);


	int total_indexFW = 0;
	int total_indexRC = 0;
	int total_Toasted = 0;
        for(int count=0;count<numranks;count++)
		total_indexFW += all_indexFW_sizes[count];

         if (!sStrand) {
            for(int count=0;count<numranks;count++)
                total_indexRC += all_indexRC_sizes[count];
	   }

	   for(int count=0;count<numranks;count++)
		total_Toasted += all_Toasted_sizes[count];


	int* FW_displs;
	int* RC_displs;
	int* Toasted_displs;

		Toasted_displs = (int*)malloc(sizeof(int)*numranks);
		Toasted_displs[0]=0;
		for(int count=1;count<numranks;count++)
		   Toasted_displs[count] = Toasted_displs[count-1] + all_Toasted_sizes[count-1]; 


		FW_displs = (int*)malloc(sizeof(int)*numranks);

		FW_displs[0] = 0;
		for(int count=1;count<numranks;count++) 
		   FW_displs[count] =  FW_displs[count-1] + all_indexFW_sizes[count-1];

		all_IndexFW_left  = (int*)malloc(sizeof(int)*total_indexFW);
		all_IndexFW_right = (int*)malloc(sizeof(int)*total_indexFW);
		all_Index_Toasted = (int*)malloc(sizeof(int)*total_Toasted); 

		MPI_Allgatherv(leftFW, numFW,MPI_INT,all_IndexFW_left ,all_indexFW_sizes, FW_displs, MPI_INT,MPI_COMM_WORLD);
		MPI_Allgatherv(rightFW,numFW,MPI_INT,all_IndexFW_right,all_indexFW_sizes, FW_displs, MPI_INT,MPI_COMM_WORLD);
		MPI_Allgatherv(Ind_Toasted,numToasted,MPI_INT,all_Index_Toasted,all_Toasted_sizes, Toasted_displs, MPI_INT,MPI_COMM_WORLD);

		for(int count=0;count<total_indexFW;count++)
                        add_reciprocal_iworm_link(weld_reinforced_iworm_clusters, all_IndexFW_left[count], all_IndexFW_right[count] );

		for(int count=0; count<total_Toasted; count++) toasted[ all_Index_Toasted[count] ] = true;	

		free(FW_displs);
		free(Toasted_displs);
		free(all_indexFW_sizes);
		free(all_Toasted_sizes);
		free(all_IndexFW_left);
		free(all_IndexFW_right);		
		free(all_Index_Toasted);
		free(leftFW);
		free(rightFW);
		free(Ind_Toasted);

		if (!sStrand) {
		  RC_displs = (int*)malloc(sizeof(int)*numranks);
		  RC_displs[0] = 0;
		  for(int count=1;count<numranks;count++)
			RC_displs[count] =  RC_displs[count-1] + all_indexRC_sizes[count-1];

		
		  all_IndexRC_left  = (int*)malloc(sizeof(int)*total_indexRC);
                  all_IndexRC_right = (int*)malloc(sizeof(int)*total_indexRC);	

		  MPI_Allgatherv(leftRC, numRC,MPI_INT,all_IndexRC_left ,all_indexRC_sizes, RC_displs, MPI_INT,MPI_COMM_WORLD);
                  MPI_Allgatherv(rightRC,numRC,MPI_INT,all_IndexRC_right,all_indexRC_sizes, RC_displs, MPI_INT,MPI_COMM_WORLD);


		  for(int count=0;count<total_indexRC;count++)
                        add_reciprocal_iworm_link(weld_reinforced_iworm_clusters, all_IndexRC_left[count], all_IndexRC_right[count] );	

		  free(RC_displs);
		  free(all_indexRC_sizes);
		  free(all_IndexRC_left);
		  free(all_IndexRC_right);
		  free(leftRC);
		  free(rightRC);

		}

	   cerr<< endl << endl << " Overhead of MPI steps  "<<timer_stop()<<" (in milliseconds)"<<endl;	


	   svec<Pool> clustered_pools = bubble_up_cluster_growth(weld_reinforced_iworm_clusters, toasted);

	    if (DEBUG) {
        	cerr << "Final pool description: " << endl;
        	describe_poolings(clustered_pools);
    		}

    		add_unclustered_iworm_contigs(clustered_pools, dna);

	MPI_Barrier(MPI_COMM_WORLD);

	    int component_count = 0;

    	  for (i=0; i<clustered_pools.isize(); i++) {
        	Pool & p = clustered_pools[i];

        	if (p.size() == 0) { continue; }

        	p.sortvec();

        	int sum_iworm_length = 0;
        	for (unsigned int j = 0; j < p.size(); j++) {
            		int z = p.get(j);
            		sum_iworm_length += dna[z].isize();
        	}

        	if (sum_iworm_length < MIN_CONTIG_LENGTH) {
            		continue;
        	}

        	stringstream pool_info;
	if(rank==0) {
        	pool_info << "#POOL_INFO\t" << component_count << ":" << "\t";

		cout << "COMPONENT " << component_count << "\t" << p.size() << endl;
        	for (unsigned int j = 0; j < p.size(); j++) {
            		int z = p.get(j);

            		pool_info << z << " ";
            		cout << ">Component_" << component_count << " " << p.size() << " " << z << " [iworm" << dna.Name(z) << "]" << endl;
            		PrintSeq(dna[z]);
        	}
        	pool_info << endl;
          
		cout << pool_info.str();

        	cout << "END" << endl;
	   }
        	component_count++;

    	  }


#else   

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
                        cerr << "Phase2: kmer match between " << dna.Name(i) << " and " << dna.Name(c)
                             << " at F positions: " << j << " and " << start << endl;
                    } 
                
                
                    //cout << "TEST_WELD: " << I << " " << J << " ";
                    string welding_kmer;
                    int welding_kmer_read_count;
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
               
                    } else if (!IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
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
                             << " with " << welding_kmer << " found in " << welding_kmer_read_count << " reads" << endl;
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
                    int welding_kmer_read_count;
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

                    } else if (!IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
                        // cerr << "Rejecting rc merge between " << dna.Name(i);
                        // cerr << " and " << dna.Name(c) << endl;
                        continue;
                    } else if (!IsGoodCoverage(coverage, coverage_other, min_iso_ratio)) {
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
                             << " with " << welding_kmer << " found in " << welding_kmer_read_count << " reads" << endl;
                    }
                
                    //cout << "Mapped sequence " << dna.NameClean(c) << " to pool " << mapped[c] << " -" << endl;
                }
            
            }
 
	}

        end_time = omp_get_wtime();
	cerr << endl << endl << "...done Phase 2. (" << end_time - start_time << " seconds)" << endl;

        cerr<<"\nPhase 2 took "<<timer_stop()<< "(in milliseconds)"<<endl;
        cerr << endl; // end of progress monitoring.
	
    
    // sl_clustered_pools = sl_cluster_pools(weld_reinforced_iworm_clusters, toasted);
    
    svec<Pool> clustered_pools = bubble_up_cluster_growth(weld_reinforced_iworm_clusters, toasted);
    
    //-----------------------------------------------------------------------------------
    // Generate final output


    if (DEBUG) {
        cerr << "Final pool description: " << endl;
        describe_poolings(clustered_pools);
    }
    
    
    add_unclustered_iworm_contigs(clustered_pools, dna);
    
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

#endif

    }
    print_cpu_memory();
#ifdef MPI_ENABLED
    MPI_Finalize();
#endif

   
   return 0;
    
}
  
