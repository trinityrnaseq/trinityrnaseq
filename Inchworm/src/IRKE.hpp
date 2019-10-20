#pragma once

#include <map>
#include <string>

#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "string_util.hpp"

using namespace std;

typedef pair<vector<kmer_int_type_t>,int> Path_n_count_pair;

struct my_sorter {
    
    bool operator() (Path_n_count_pair& a, Path_n_count_pair& b) {
        
        if (a.second < b.second) {
            return(true);
        }
        else {
            return(false);
        }
    }
}; 


class IRKE {
    
public:
    
    IRKE();
    IRKE(unsigned int kmer_length, 
         unsigned int max_recursion = 25,
         float min_seed_entropy = 1.5,
         unsigned int min_seed_coverage = 2,
         float min_any_entropy = 0.0,
         bool pacman = false,
         bool crawl = false,
         unsigned int crawl_length = 1,
         bool double_stranded = false
         );


    bool TWO_PHASE; // draft inchworm contig built, new seed selected from it, and then proper inchworm contig built.
    
    void build_graph(const string& reads_fasta_filename, bool reassembleIworm = false, bool useKmers = false);
    
    unsigned long get_graph_size();
    
    void set_prune_singleton_read_interval (unsigned long interval);
    
    void prune_kmers_min_count(unsigned int min_count);
    
    void prune_kmers_min_entropy(float min_entropy);
    
    bool prune_kmer_extensions(float min_ratio_non_error);
    
    bool prune_some_kmers(unsigned int min_count, float min_entropy, bool prune_error_kmers, float min_ratio_non_error);

    bool prune_branched_kmers();
    
    //void report_connected_components (float min_connectivity, unsigned int min_component_size, bool summarize_size_only_flag = false);
    
    string reconstruct_path_sequence(vector<kmer_int_type_t>& path, vector<unsigned int>& cov_counter);
    
    bool sequence_path_exists(string& sequence, unsigned int min_coverage, float min_entropy, float min_connectivity, vector<unsigned int>& cov_counter);
    
    void monitor_progress(unsigned int monitor_setting);
    
    void compute_sequence_assemblies(float min_connectivity = 0.2, 
                                     unsigned int min_assembly_length = 100, 
                                     unsigned int min_assembly_coverage = 2,
                                     bool PARALLEL_IWORM = false,
                                     bool write_coverage = false, 
                                     string coverage_output_filename = "assembly_coverage.out");
    
    
    void describe_kmers(); // just describe each kmer in the graph, print to stdout.
    
    string thread_sequence_through_graph(string& sequence); // describe kmers in this sequence.
    

    // utilities for setting and unsetting the sorted kmer list targeted as seeds for inchworm assembly
    void populate_sorted_kmers_list();
    void clear_sorted_kmers_list();
    
    
    
private:
    
    KmerCounter kcounter;
    
    unsigned int MONITOR;
    unsigned int INCHWORM_ASSEMBLY_COUNTER;
    float MIN_SEED_ENTROPY;
    unsigned int MIN_SEED_COVERAGE;
    float MIN_ANY_ENTROPY;
    bool PACMAN;
    bool CRAWL;
    unsigned int CRAWL_LENGTH;
    unsigned int MAX_RECURSION;
    bool DOUBLE_STRANDED_MODE;
    unsigned long PRUNE_SINGLETON_READ_INTERVAL;

    

    bool got_sorted_kmers_flag;
    //vector<Kmer_counter_map_iterator> sorted_kmers;
    vector<Kmer_Occurence_Pair> sorted_kmers;
    
    void populate_Kmers_from_kmers(const string& kmers_filename);
    void populate_Kmers_from_fasta(const string& fasta_filename, bool reassembleIworm = false);
    
    string reconstruct_path_sequence(KmerCounter& kcounter, vector<kmer_int_type_t>& path, vector<unsigned int>& cov_counter);
    
    
    Path_n_count_pair inchworm(KmerCounter& kcounter, char direction, kmer_int_type_t kmer,
                               Kmer_visitor& visitor, float min_connectivity);  // direction = [FR]
    
    Path_n_count_pair inchworm_step(KmerCounter& kcounter, char direction, Kmer_Occurence_Pair kmer, Kmer_visitor& visitor,
                                    Kmer_visitor& eliminator, unsigned int inchworm_round, unsigned int depth, 
                                    float min_connectivity, unsigned int max_recurse);  // direction = [FR]
    
    void compute_sequence_assemblies(KmerCounter& kcounter, float min_connectivity, 
                                     unsigned int min_assembly_length, unsigned int min_assembly_coverage,
                                     bool PARALLEL_IWORM, bool write_coverage, string coverage_output_filename);
    
    void traverse_path(KmerCounter& kcounter, Kmer_Occurence_Pair seed_kmer, Kmer_visitor& visitor, Kmer_visitor& place_holder,
                       float min_connectivity_ratio, unsigned int depth);
    
    // helper function prototypes
    vector<kmer_int_type_t> _join_forward_n_reverse_paths(vector<kmer_int_type_t>& reverse_path, 
                                                          kmer_int_type_t seed_kmer_val, vector<kmer_int_type_t>& forward_path);
    
    bool exceeds_min_connectivity (KmerCounter& kcounter, Kmer_Occurence_Pair kmerA, Kmer_Occurence_Pair kmerB, float min_connectivity);
    bool exceeds_min_connectivity (KmerCounter& kcounter, string kmerA, string kmerB, float min_connectivity);
    
    void describe_kmers(KmerCounter& kcounter); 
        
    vector<kmer_int_type_t> build_inchworm_contig_from_seed(kmer_int_type_t kmer, KmerCounter& kcounter, 
                                                            float min_connectivity, unsigned int& total_counts,
                                                            bool PARALLEL_IWORM);
    
    bool is_good_seed_kmer(kmer_int_type_t kmer, unsigned int kmer_count, unsigned int kmer_length,
                           float min_connectivity);
    
    kmer_int_type_t extract_best_seed(vector<kmer_int_type_t>& kmer_vec, KmerCounter& kcounter, float min_connectivity);
        
    
};

