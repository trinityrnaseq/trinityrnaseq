#ifndef __KMERCOUNTER__
#define __KMERCOUNTER__

#include <vector>
#include <map>
#include <string>
#include <set>


#include "sequenceUtil.hpp"


using namespace std;



struct eqstr {

  bool operator()(const kmer_int_type_t& kmer_val_a, const kmer_int_type_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};


struct hashme {
  
  kmer_int_type_t operator() (const kmer_int_type_t& kmer_val) const {

	return(kmer_val);
  }
};


#ifdef __GOOGLE__

// #warning "******** using GOOGLE SPARSEHASH for Kmer graph *********"

#include <google/sparse_hash_map>
using google::sparse_hash_map; 

typedef sparse_hash_map<kmer_int_type_t, unsigned int, hashme, eqstr> Kmer_counter_map;
typedef sparse_hash_map<kmer_int_type_t, unsigned int, hashme, eqstr>::iterator Kmer_counter_map_iterator;
typedef sparse_hash_map<kmer_int_type_t, unsigned int, hashme, eqstr>::const_iterator Kmer_counter_map_const_iterator;

#elif defined(__SUNPRO_CC) // Solaris Studio compiler

#include <hash_map>
typedef std::hash_map<kmer_int_type_t,unsigned int, hashme, eqstr> Kmer_counter_map;
typedef std::hash_map<kmer_int_type_t,unsigned int, hashme, eqstr>::iterator Kmer_counter_map_iterator;
typedef std::hash_map<kmer_int_type_t,unsigned int, hashme, eqstr>::const_iterator Kmer_counter_map_const_iterator;

#else
#include <ext/hash_map>

typedef __gnu_cxx::hash_map<kmer_int_type_t,unsigned int, hashme, eqstr> Kmer_counter_map;
typedef __gnu_cxx::hash_map<kmer_int_type_t,unsigned int, hashme, eqstr>::iterator Kmer_counter_map_iterator;
typedef __gnu_cxx::hash_map<kmer_int_type_t,unsigned int, hashme, eqstr>::const_iterator Kmer_counter_map_const_iterator;

#endif

typedef pair<kmer_int_type_t,unsigned int> Kmer_Occurence_Pair;


class KmerCounter {
  
public:
  
    KmerCounter() {};
    KmerCounter(unsigned int kmer_length, bool is_ds=false);
    
    unsigned int get_kmer_length();
    unsigned long size();
    
    void add_sequence(string& sequence, unsigned int cov=1);
    
    bool add_kmer (kmer_int_type_t, unsigned int count);
    bool add_kmer (string kmer, unsigned int count);
    
    void describe_kmers();
    void dump_kmers_to_file(string& outfilename);
    
    //vector<Kmer_counter_map_iterator> get_kmers_sort_descending_counts();
    vector<Kmer_Occurence_Pair> get_kmers_sort_descending_counts();
    
    Kmer_counter_map_iterator find_kmer(kmer_int_type_t kmer_val);
    
    bool kmer_exists(string kmer);
    bool kmer_exists(kmer_int_type_t kmer_val);
    
    unsigned int get_kmer_count(string kmer);
    unsigned int get_kmer_count(kmer_int_type_t kmer_val);
    
    string describe_kmer(string& kmer);
    
    string get_kmer_string(kmer_int_type_t kmer_val);
    kmer_int_type_t get_kmer_intval(string kmer);
    
    
    bool prune_kmer(string kmer); // remove kmer from map
    bool prune_kmer(kmer_int_type_t kmer_val);
    bool prune_some_kmers(unsigned int min_count, float min_entropy, bool prune_error_kmers, float min_ratio_non_error);

    bool prune_branched_kmers();
    
    void prune_kmers_min_count(unsigned int count);
    void prune_kmers_min_entropy(float min_entropy);
    
    bool clear_kmer(kmer_int_type_t kmer_val);
    
    // methods return kmers sorted descending by count.
    //vector<string> get_reverse_kmer_candidates(string& kmer);
    vector<Kmer_Occurence_Pair> get_reverse_kmer_candidates(kmer_int_type_t seed_kmer);
    //vector<string> get_forward_kmer_candidates(string& kmer);
    vector<Kmer_Occurence_Pair> get_forward_kmer_candidates(kmer_int_type_t seed_kmer);
    
    // methods return kmers unsorted, in order G,A,T,C
    //vector<string> get_reverse_kmer_candidates_unsorted(string& kmer);
    vector<Kmer_Occurence_Pair> get_reverse_kmer_candidates_unsorted(kmer_int_type_t seed_kmer, bool getZeros);
    //vector<string> get_forward_kmer_candidates_unsorted(string& kmer);
    vector<Kmer_Occurence_Pair> get_forward_kmer_candidates_unsorted(kmer_int_type_t seed_kmer, bool getZeros);
    
    // get the simple list of the 4 possible kmer extensions.
    kmer_int_type_t* get_forward_kmer_candidates_noLookup(kmer_int_type_t seed_kmer, kmer_int_type_t forward_kmer_array_size_4 [4]);    
    kmer_int_type_t* get_reverse_kmer_candidates_noLookup(kmer_int_type_t seed_kmer, kmer_int_type_t reverse_kmer_array_size_4 [4]);

    
    bool prune_kmer_extensions( float min_ratio_non_error);
    
    const Kmer_counter_map& get_kmer_counter_map() const;
    
  
  
private:
  
  unsigned int _kmer_length;
  
  Kmer_counter_map _kmer_counter;

  bool _DS_MODE;

};



class Kmer_visitor {
public:
	Kmer_visitor(unsigned int kmer_length, bool is_ds);
	void add (kmer_int_type_t kmer);
	bool exists (kmer_int_type_t kmer);
	void erase (kmer_int_type_t kmer);
	void clear();
	unsigned int size();

private:
	unsigned int _kmer_length;
	bool _DS_MODE;
	set<kmer_int_type_t> _set;
};



class Sort_kmer_by_count_desc {
  
public:
  
  Sort_kmer_by_count_desc(KmerCounter *kcounter);
  
  bool operator() (const Kmer_counter_map_iterator& i, const Kmer_counter_map_iterator& j);

  bool operator() (const Kmer_Occurence_Pair& i, const Kmer_Occurence_Pair& j);

  bool operator() (const kmer_int_type_t& val_i, const kmer_int_type_t& val_j);
  
  bool operator() (const string& i, const string& j);

private:
  
  KmerCounter *kcounter;

};



#endif
