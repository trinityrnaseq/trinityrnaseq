#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include "IRKE.hpp"
#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "MaskedKmerCounter.hpp"
#include "stacktrace.hpp"
#include "irke_common.hpp"
#include "argProcessor.hpp"

unsigned int IRKE_COMMON::MONITOR = 0;
size_t KMER_SIZE = 25;
int MAX_THREADS = 6;

// various devel params
bool IRKE_COMMON::__DEVEL_no_kmer_sort = false;

void populate_kmer_counter_from_kmers (MaskedKmerCounter& kcounter, string& kmers_fasta_file);
string usage();


string stomp_snps(string& sequence, MaskedKmerCounter& kcounter, unsigned long& snp_stomp_counter);

int runMe(int argc, char* argv[]);

int main(int argc, char* argv[]) {

    try {
        return(runMe(argc, argv));
    }
    catch (string errmsg) {
        
        cerr << "ERROR encountered: " << errmsg;
        exit(1);
        
    }

    return(0);
}

int runMe(int argc, char* argv[]) {

    ArgProcessor args(argc, argv);
    if(args.isArgSet("--help") ||
       (!(args.isArgSet("--reads")
          &&
          args.isArgSet("--kmers") )
          ) ) {
        cerr << usage() << endl << endl;
        exit(1);
    }
    
    string reads_fasta_file = args.getStringVal("--reads");
    
    bool is_DS = (! args.isArgSet("--SS"));
    if(args.isArgSet("--kmer_size")) {
        KMER_SIZE = args.getIntVal("--kmer_size");
        if(KMER_SIZE < 20) {
            cerr << "Error, min kmer size is 20";
            exit(2);
        }
        if (KMER_SIZE % 2 != 1) {
            cerr << "Error, kmer size: " << KMER_SIZE << " is not an odd number.  Must be odd." << endl;
            exit(3);
        }
    }
    if(args.isArgSet("--monitor")) {
        IRKE_COMMON::MONITOR = args.getIntVal("--monitor");
     }
    if (args.isArgSet("--num_threads")) {
        int num_threads = args.getIntVal("--num_threads");
        if (num_threads < MAX_THREADS) {
            omp_set_num_threads(num_threads);
        }
        else {
            // set to max
            omp_set_num_threads(MAX_THREADS);
        }
    }
    
    if(omp_get_max_threads() > MAX_THREADS) {
        omp_set_num_threads(MAX_THREADS);
    }
    MaskedKmerCounter kcounter (KMER_SIZE, is_DS);
    
    string kmers_fasta_file = args.getStringVal("--kmers");
    populate_kmer_counter_from_kmers(kcounter, kmers_fasta_file);

    //kcounter.describe_kmers();
    //exit(2);
        
    int start_time = time(NULL);

    Fasta_reader fasta_reader(reads_fasta_file);

    unsigned long snp_stomp_counter = 0;
    
    while (true) {

        if (! fasta_reader.hasNext())
            break;

        Fasta_entry fe = fasta_reader.getNext();
        string sequence = fe.get_sequence();
        if(sequence == "")
            continue;

        string header = fe.get_header();
        
        
        

        string newseq = stomp_snps(sequence, kcounter, snp_stomp_counter);


        cout << ">" << header << endl << newseq << endl;

                
    }
    
    int end_time = time(NULL);

    cerr << "Number of snps stomped: " << snp_stomp_counter << endl;
    
    cerr << "SNP_STOMPING_TIME: " << (end_time - start_time) << " seconds." << endl;

    
    return(0);
}

void populate_kmer_counter_from_kmers(MaskedKmerCounter& kcounter, string& kmers_fasta_file) {

    unsigned long start, end;
    cerr << "-reading Kmer occurrences..." << endl;
    start = time(NULL);
    Fasta_reader fasta_reader(kmers_fasta_file);

    unsigned int record_counter = 0;
    while (true) {
        Fasta_entry fe = fasta_reader.getNext();
        if(fe.get_sequence() == "") break;
        if(IRKE_COMMON::MONITOR) {
            if(record_counter % 100000 == 0)
                cerr << "\r [" << record_counter/1000000 << "M] Kmers parsed.     ";
        }
        string seq = fe.get_sequence();
        if(seq.length() != KMER_SIZE) {
            cerr << "ERROR: kmer " << seq << " is not of length: " << KMER_SIZE << endl;
            continue;
        }
        record_counter += 1;
        kmer_int_type_t kmer_val = kcounter.get_kmer_intval(seq);
        unsigned int count = atoi(fe.get_header().c_str());
        kcounter.set_kmer_occurrence_pair(kmer_val, count);
    }
    
    end = time(NULL);
    cerr << endl << " done parsing " << record_counter << " Kmers, " << kcounter.size() << " added, taking " << (end-start) << " seconds." << endl;
    return;
}


string stomp_snps(string& sequence, MaskedKmerCounter& kcounter, unsigned long& snp_stomp_counter) {
    if(IRKE_COMMON::MONITOR) {
        cerr << "processing sequence: " << sequence << endl;
    }
    if (sequence.length() < KMER_SIZE) {
        return(sequence);
    }
    

    string sequence_copy = sequence;

    unsigned int central_kmer_pos = KMER_SIZE / 2;
    
    for (size_t i = 0; i <= sequence.length() - KMER_SIZE; i++) {
        // cerr << "i: " << i << ", <= " << sequence.length() - KMER_SIZE << endl;
        string kmer = sequence.substr(i, KMER_SIZE);
                        
        if(contains_non_gatc(kmer)) {
            continue;
        }
        unsigned int string_pos = i + central_kmer_pos;

        kmer_int_type_t kmer_val = kmer_to_intval(kmer);
        kmer_int_type_t ds_kmer_int_val = get_DS_kmer_val(kmer_val, KMER_SIZE);


        kmer_int_type_t masked_kmer_int_val = kcounter.get_masked_kmer_intval(kmer_val);
        kmer_int_type_t masked_ds_kmer_int_val = kcounter.get_masked_kmer_intval(ds_kmer_int_val);

        kmer_int_type_t search_kmer_val = (kcounter.is_DS_mode()) ? masked_ds_kmer_int_val : masked_kmer_int_val;
        string search_kmer = decode_kmer_from_intval(search_kmer_val, KMER_SIZE);
        
        Kmer_Occurrence_Pair kmer_occur_pair = kcounter.get_kmer_occurrence_pair(search_kmer_val);

        kmer_int_type_t best_unmasked_kmer_val = kmer_occur_pair.first;

        unsigned int count = kmer_occur_pair.second;

        string stored_kmer = decode_kmer_from_intval(best_unmasked_kmer_val, KMER_SIZE);

        if (IRKE_COMMON::MONITOR) {
            cerr << search_kmer << "\t" << stored_kmer << "\t" << count << endl;
        }
        char chosen_char = stored_kmer[central_kmer_pos];
        if (kcounter.is_DS_mode()) {
            // see if stored kmer matches the forward or reverse of the read kmer
            kmer_int_type_t masked_stored_val = kcounter.get_masked_kmer_intval(best_unmasked_kmer_val);
            if (masked_stored_val == masked_kmer_int_val) {
                // forward orientation
            }
            else if (masked_stored_val == masked_ds_kmer_int_val) {
                // reverse orientation
                string schar = "";
                schar += chosen_char;
                schar = revcomp(schar);
                chosen_char = schar[0];
            }
        }
        
        if (chosen_char != sequence[string_pos]) {

            if (IRKE_COMMON::MONITOR) {
                stringstream ss;
                ss << "snp [" << string_pos << "]: " << sequence[string_pos] << " -> " << chosen_char << endl;
                cerr << ss.str();
            }

            // update sequence w/ stomped snp
            sequence_copy[string_pos] = chosen_char;

            snp_stomp_counter += 1;
            
        }
    }
    
    return(sequence_copy);
}


string usage() {
    stringstream usage_info;
    usage_info << endl << endl;
    usage_info << "Usage: " << endl
               << "  --reads  <str>             " << ":fasta file containing target reads for kmer coverage stats" << endl
               << "  --kmers  <str>             " << ":fasta file containing kmers" << endl
               << endl 
               << "* optional:" << endl
               << "  --kmer_size <int>          " << ":kmer size, must be odd value.  default = 25" << endl
               << "  --SS                       " << ":strand specific mode.  (default: double-stranded RNA-Seq mode (not strand-specific))" << endl
               << "  --capture_coverage_info    " << ":writes coverage info file." << endl
               << "  --monitor <int>            " << ":verbose level for debugging" << endl
               << "  --num_threads <int>        " << ":number of threads" << endl
               << endl;
    return(usage_info.str());
}

