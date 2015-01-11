#include "Fasta_reader.hpp"
#include "Fasta_entry.hpp"
#include "KmerCounter.hpp"
#include <algorithm>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <iostream>
#include "argProcessor.hpp"
#include "sequenceUtil.hpp"
#include "irke_common.hpp"

void parse_target_kmers(string& filename, KmerCounter& kcounter);
bool contains_kmer(string& sequence, KmerCounter& kcounter, unsigned int kmer_length, 
                   int max_kmer_count, float min_kmer_entropy, int min_percent_read_containing_kmers);
int run(int argc, char* argv[]);
void print_kmer_count_histogram(KmerCounter& kcounter);


// various devel params
bool IRKE_COMMON::__DEVEL_no_kmer_sort = false;

int main (int argc, char* argv[]) {
    
    try {
        return(run(argc, argv));
    }
    catch (string err) {
        cerr << err;
        exit(1);
    }
    
}

int run (int argc, char* argv[]) {
    
    ArgProcessor args(argc, argv);
    
    
    string targets_fa_filename;
    string reads_fa_filename;
    
    unsigned int kmer_length = 25;
    bool print_kmer_hist = false;
    bool DS_mode = false;
    int max_kmer_count = -1;
    float min_kmer_entropy = 1.5f;
    int min_percent_read_containing_kmers = -1;
    
    stringstream usage;
    
    usage << "usage: " << argv[0] << " --target seq_or_kmers.fasta  --reads reads.fasta [options]" 
          << endl << endl
          << "--print_kmer_hist                          prints histogram of kmer abundances" << endl
          << "--DS                                       treats target.fasta as double-stranded" << endl
          << "-K <int>                                   Kmer size (default: 25) " << endl
          << "--max_kmer_count <int>                     maximum occurrence to use kmer for read recruitment" << endl
          << "--min_kmer_entropy <float>                 minimum kmer entropy to use for recruitment (default: 1.5)" << endl
          << "--min_percent_read_containing_kmers <int>  minimum percent of the read that is covered by highly-occurring kmers" << endl
          << endl << endl;
    
    
    // process options:
    if (! (args.isArgSet("--target") 
           && 
           (args.isArgSet("--reads") || args.isArgSet("--print_kmer_hist") )
           )
        ) {
        
        cerr << usage.str();
        return(1);
    }
    else {
        
        targets_fa_filename = args.getStringVal("--target");
        reads_fa_filename = args.getStringVal("--reads");
        
    }
    
    if (args.isArgSet("--print_kmer_hist")) {
        print_kmer_hist = true;
    }
    
    if (args.isArgSet("--DS")) {
        DS_mode = true;
    }
    
    if (args.isArgSet("-K")) {
        kmer_length = args.getIntVal("-K");
    }
    
    if (args.isArgSet("--max_kmer_count")) {
        max_kmer_count = args.getIntVal("--max_kmer_count");
    }
    
    if (args.isArgSet("--min_kmer_entropy")) {
        min_kmer_entropy = args.getFloatVal("--min_kmer_entropy");
    }
    
    if (args.isArgSet("--min_percent_read_containing_kmers")) {
        min_percent_read_containing_kmers = args.getIntVal("--min_percent_read_containing_kmers");
    }
    
    
    KmerCounter kcounter(kmer_length, DS_mode);
    
    parse_target_kmers(targets_fa_filename, kcounter);
    
    if (print_kmer_hist) {
        print_kmer_count_histogram(kcounter);
        return(0);
    }
    
    unsigned int counter = 0;
    
    Fasta_reader fr(reads_fa_filename);
    while (fr.hasNext()) {
        
        counter++;
        if (counter % 1000 == 0) {
            cerr << "\r[" << counter << "]  reads parsed      ";
        }
        
        
        Fasta_entry fe = fr.getNext();
        string seq = fe.get_sequence();
        
        if (contains_kmer(seq, kcounter, kmer_length, max_kmer_count, min_kmer_entropy, min_percent_read_containing_kmers)) {
            cout << ">" << fe.get_header() << endl << seq << endl;
        }
    }
    
    
    return(0);
}


void parse_target_kmers(string& filename, KmerCounter& kcounter) {
    
    Fasta_reader fr(filename);
    
    while (fr.hasNext()) {
        Fasta_entry fe = fr.getNext();
        string seq = fe.get_sequence();
        
        // cerr << "parsing sequence: " << seq << endl;
        
        kcounter.add_sequence(seq);
    }
    
    return;
}

bool contains_kmer(string& sequence, KmerCounter& kcounter, unsigned int kmer_length, int max_kmer_count, 
                   float min_kmer_entropy, int min_percent_read_containing_kmers) {
    
    
    int num_kmers_contained = 0;
    for (unsigned int i = 0; i <= sequence.length() - kmer_length; i++) {
        
        string kmer = sequence.substr(i, kmer_length);
        
        if ((! contains_non_gatc(kmer)) && kcounter.kmer_exists(kmer)) {
            
            // check kmer abundance
            int count = kcounter.get_kmer_count(kmer);
            if (max_kmer_count > 0 && count > max_kmer_count) {
                continue;
            }
            
            // check kmer sequence complexity
            float entropy = compute_entropy(kmer);
            if (entropy < min_kmer_entropy) {
                continue;
            }
            
            if (min_percent_read_containing_kmers > 0) {
                num_kmers_contained++;
            }
            else {
                return(true); // one kmer is enough to suffice
            }
        }
        
        
    }

    if (min_percent_read_containing_kmers > 0) {
        float percent_contained_by_kmers = (float) num_kmers_contained / (sequence.length() - kmer_length) * 100;
   
        if (percent_contained_by_kmers >= min_percent_read_containing_kmers) {
            return(true);
        }
    }
    
    
    return(false); // kmer doesn't exist
}


void print_kmer_count_histogram(KmerCounter& kcounter) {
    
    cerr << "generating kmer histogram" << endl;
    
    map<unsigned int,unsigned int> count_counter;
    
    const Kmer_counter_map& kmer_map = kcounter.get_kmer_counter_map();
    
    Kmer_counter_map_const_iterator it;
    
    for (it = kmer_map.begin(); it != kmer_map.end(); it++) {
        
        unsigned int count = it->second;
        if (count_counter.find(count) == count_counter.end()) {
            count_counter[count] = 1;
        }
        else {
            count_counter[count] = count_counter.find(count)->second + 1;
        }
                
    }
    
    cerr << "done counting kmers." << endl;
    
    map<unsigned int,unsigned int>::iterator it2;
    vector<unsigned int> counts;
    
    
    for (it2 = count_counter.begin(); it2 != count_counter.end(); it2++) {
        counts.push_back(it2->first);
    }
    
    
    cerr << "Sorting kmer counts .... ";
    sort(counts.begin(), counts.end());
    cerr << "done. " << endl;
    
    vector<unsigned int>::iterator it3;
    for (it3 = counts.begin(); it3 != counts.end(); it3++) {
        
        unsigned int count = *it3;
        unsigned int kmer_counts = count_counter[count];
        
        cout << count << "\t" << kmer_counts << endl;
        //unsigned int 
        
    }
    
    cout << endl;
    
    return;
}
