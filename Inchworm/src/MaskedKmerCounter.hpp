#ifndef __MASKEDKMERCOUNTER__
#define __MASKEDKMERCOUNTER__

#include <vector>
#include <map>
#include <string>
#include <set>


#include "sequenceUtil.hpp"


using namespace std;



struct MaskedKmerCounter__eqstr {

  bool operator()(const kmer_int_type_t& kmer_val_a, const kmer_int_type_t& kmer_val_b) const {
	return( kmer_val_a == kmer_val_b);
  }
};


struct MaskedKmerCounter__hashme {
  
  kmer_int_type_t operator() (const kmer_int_type_t& kmer_val) const {

	return(kmer_val);
  }
};

#include <ext/hash_map>

typedef pair<kmer_int_type_t,unsigned int> Kmer_Occurrence_Pair;

typedef __gnu_cxx::hash_map<kmer_int_type_t, Kmer_Occurrence_Pair, MaskedKmerCounter__hashme, MaskedKmerCounter__eqstr> Masked_kmer_map;
typedef __gnu_cxx::hash_map<kmer_int_type_t, Kmer_Occurrence_Pair, MaskedKmerCounter__hashme, MaskedKmerCounter__eqstr>::iterator Masked_kmer_map_iterator;
typedef __gnu_cxx::hash_map<kmer_int_type_t, Kmer_Occurrence_Pair, MaskedKmerCounter__hashme, MaskedKmerCounter__eqstr>::const_iterator Masked_kmer_map_const_iterator;


class MaskedKmerCounter {
  
public:
  
    MaskedKmerCounter() {};
    MaskedKmerCounter(unsigned int kmer_length, bool is_ds=false);
    
    unsigned int get_kmer_length();
    unsigned long size();

    void describe_kmers();

    void dump_kmers_to_file(string& outfilename);

    Masked_kmer_map_iterator find_kmer(kmer_int_type_t kmer_val);

    bool kmer_exists(string kmer);
    bool kmer_exists(kmer_int_type_t kmer_val);
    
    unsigned int get_kmer_count(string kmer);
    unsigned int get_kmer_count(kmer_int_type_t kmer_val);

    string get_kmer_string(kmer_int_type_t kmer_val);    

    Kmer_Occurrence_Pair get_kmer_occurrence_pair(kmer_int_type_t kmer_val);
    Kmer_Occurrence_Pair get_kmer_occurrence_pair(string kmer);
    
    kmer_int_type_t get_kmer_intval(string kmer);
    
    kmer_int_type_t get_masked_kmer_intval(kmer_int_type_t kmer_val);

    bool set_kmer_occurrence_pair(kmer_int_type_t kmer_val, unsigned int count);
    bool set_kmer_occurrence_pair(string kmer, unsigned int count);
    
    string describe_kmer(string& kmer);
    string describe_kmer(kmer_int_type_t kmer);
    
    const Masked_kmer_map& get_masked_kmer_map() const;

    bool is_DS_mode();
  
  
private:
  
    unsigned int _kmer_length;
    unsigned int _masked_pos;
    kmer_int_type_t  _masked_bit_fields;
    
    Masked_kmer_map _masked_kmer_map;
    
    bool _DS_MODE;
    
};

#endif
