#pragma once

#include <map>
#include "analysis/sequenceUtil.h"




class DeBruijnKmer {
    
public:

    DeBruijnKmer(kmer_int_type_t k);
    DeBruijnKmer(const DeBruijnKmer& dk);
    
    unsigned int getID() const;
    kmer_int_type_t get_kmer_int_val() const;
    
    vector<kmer_int_type_t> get_prev_kmers(unsigned int kmer_length);
    vector<kmer_int_type_t> get_next_kmers(unsigned int kmer_length);
    
    void add_prev_kmer(kmer_int_type_t k, unsigned int kmer_length);
    void add_next_kmer(kmer_int_type_t k, unsigned int kmer_length);
    
    
    //private:
    
    static unsigned int id_counter;
    unsigned int id;
    kmer_int_type_t _kmer ;
    
    char _prev; // bit array GATC indicating prev kmers
    char _next; // ditto for next kmers
    
    static const char _G_mask; // = 8;
    static const char _A_mask; // = 4;
    static const char _T_mask; // = 2;
    static const char _C_mask; // = 1;
    

    
    
};


typedef map<kmer_int_type_t,DeBruijnKmer> DeBruijnKmerMap;


class DeBruijnGraph {

public:

    DeBruijnGraph(unsigned int kmer_length);

    void add_sequence(const string& sequence);
    
    DeBruijnKmer& get_kmer_node(kmer_int_type_t t);    
    
    string toString();
    
    string toChrysalisFormat(int component_id, bool sStrand);
    
    vector<kmer_int_type_t> get_root_kmers();
    
    bool kmerExists(kmer_int_type_t kval);
    
    vector<string> get_candidate_weldmers(kmer_int_type_t kmer_val, int weldmer_length);
    
    unsigned int get_kmer_length();
    

private:

    unsigned int _kmer_length;
    
    map<kmer_int_type_t,DeBruijnKmer> _kmer_map;
    
        
    vector<DeBruijnKmer> deconvolute_DS_mirror_graph();
    
    void recursively_construct_kmer_extensions(kmer_int_type_t seed_kmer_val, 
                                               vector<char>& kmer_extension_chars, 
                                               vector<string>& extension_kmer_strings, 
                                               char direction, // L | R  fg left or right from the seed.
                                               map<kmer_int_type_t,bool>& kmer_seen,
                                               int flank_extension_length);            
    
    
    
    string get_kmer_from_char_vector(vector<char> char_vec, char direction);
    
};
