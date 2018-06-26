#pragma once

#include <map>
#include "sequenceUtil.hpp"




class DeBruijnKmer {
    
public:

    DeBruijnKmer(kmer_int_type_t k, long long kmer_id); // kmer_id is just a unique identifier for that node, having nothing inherrently to do w/ sequence info itself.
    DeBruijnKmer(const DeBruijnKmer& dk);
    
    long long getID() const;
    kmer_int_type_t get_kmer_int_val() const;
    
    vector<kmer_int_type_t> get_prev_kmers(unsigned int kmer_length);
    vector<kmer_int_type_t> get_next_kmers(unsigned int kmer_length);
    
    void add_prev_kmer(kmer_int_type_t k, unsigned int kmer_length);
    void add_next_kmer(kmer_int_type_t k);

    string toString(int kmer_length);
    
    unsigned int increment_kmer_count(unsigned int kmer_count);
    unsigned int get_kmer_count();
    
    //private:
    
    long long id;
    kmer_int_type_t _kmer ;
    unsigned int _kmer_count;
    
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

    void add_sequence(const string& sequence, bool sStrand, unsigned int cov_val);
    
    DeBruijnKmer& get_kmer_node(kmer_int_type_t t);    
    
    string toString();
    string toDOT(bool sStrand);    
    string toChrysalisFormat(int component_id, bool sStrand);
    
    vector<kmer_int_type_t> get_root_kmers(bool sStrand);
    
    bool kmerExists(kmer_int_type_t kval);
    
    vector<string> get_candidate_weldmers(kmer_int_type_t kmer_val, int weldmer_length);
    
    unsigned int get_kmer_length();
    

private:

    unsigned int _kmer_length;
    long long _kmer_id_counter;


    map<kmer_int_type_t,DeBruijnKmer> _kmer_map;
    
        
    void recursively_construct_kmer_extensions(kmer_int_type_t seed_kmer_val, 
                                               vector<char>& kmer_extension_chars, 
                                               vector<string>& extension_kmer_strings, 
                                               char direction, // L | R  fg left or right from the seed.
                                               map<kmer_int_type_t,bool>& kmer_seen,
                                               int flank_extension_length);            
    
    
    
    string get_kmer_from_char_vector(vector<char> char_vec, char direction);
    
};
