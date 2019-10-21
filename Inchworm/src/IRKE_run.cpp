#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <omp.h>

#include "IRKE.hpp"
#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"
#include "argProcessor.hpp"
#include "irke_common.hpp"


// IRKE = Inchworm Recursive Kmer Extension  (what I called it before naming it simply 'inchworm')
//         .... the concept behind Inchworm was borrowed from Manfred Grabherr (co-author), as similarly
//         implemented in his Ananas tool.  
//                                         (bhaas)


using namespace std;


int IRKE_COMMON::MAX_THREADS_READ_PARSING = 6;  // going higher leads to decreased performance in Inchworm due to thread collisions.
int IRKE_COMMON::NUM_THREADS = -1; // use OMP_NUM_THREADS by default.
bool IRKE_COMMON::KEEP_TMP_FILES = false;

unsigned int IRKE_COMMON::MONITOR = 0; // shared among rest of code via irke_common.hpp, but lives here.

// various devel params
bool IRKE_COMMON::__DEVEL_no_kmer_sort = false;
bool IRKE_COMMON::__DEVEL_no_greedy_extend = false;
bool IRKE_COMMON::__DEVEL_no_tie_breaking = false;
bool IRKE_COMMON::__DEVEL_zero_kmer_on_use = false;
bool IRKE_COMMON::__DEVEL_rand_fracture = false;
float IRKE_COMMON::__DEVEL_rand_fracture_prob = 0; // set to probability of fracturing at any inchworm contig extension



int run_IRKE(int argc, char* argv[]);

void examine_CDS_paths(IRKE& irke, string cds_fasta_filename, unsigned int min_cov, float min_connectivity, float min_entropy,
					   bool write_coverage_flag, string coverage_filename);


void thread_sequences_through_graph(IRKE& irke, string cds_fasta_filename);

vector<string> parse_file_listing(string filename); // reads in a list of filenames from file

string usage(ArgProcessor args);

int main (int argc, char* argv[]) {
    
    try {
        return(run_IRKE(argc, argv));
    }
    
    catch (string err) {
        cerr << err << endl;
    }
    
    return(1);
    
}


int run_IRKE(int argc, char* argv[]) {

    int prog_start = time(NULL);
        
    string progname (argv[0]);
    
    // defaults:
    unsigned int MAX_RECURSION = 1;  // inchworm step size.
    
    
    bool describe_kmers = false; // if true, just dumps a report of all kmers in the graph and then exits.

    // sequence data input file options
    string reads_fasta_filename;
    string kmers_fasta_filename; // single file of k-mers (such as from Meryl)
    string kmer_files_list_filename;  // a file that contains names of other files that contain the k-mers.

    unsigned int kmer_length = 25;
    unsigned int min_kmer_count = 1;
    float MIN_CONNECTIVITY_RATIO = 0.0;
    float min_any_entropy = 0.0;
    
    bool double_stranded_flag = false;  // strand-specific by default
    
    /*************************************/
    // Report inchworm assemblies option:
    bool run_inchworm = false;
    unsigned int MIN_ASSEMBLY_LENGTH = kmer_length; // minimum length of an inchworm assembly for reporting.
    unsigned int MIN_ASSEMBLY_COVERAGE = 2; // minimum average kmer coverage for assembly to be reported.
    float MIN_SEED_ENTROPY = 1.5; // minimum entropy for a Kmer to seed in inchworm assembly construction.
    unsigned int MIN_SEED_COVERAGE = 2; // minimum kmer coverage for a seed.
    
    bool PACMAN = false;  // disables revisiting nodes previously visited within an inchworm session
    bool CRAWL = false;
    unsigned int crawl_length = 1;
    
    bool reassembleIworm = false;
    unsigned long pruneSingletonKmersReadInterval = 0;
    
    /***************************************/
    // search for CDS in the graph options:
    bool search_CDS_flag = false;
    bool thread_CDS_flag = false;
    string cds_fasta_filename;
    
    
    
    // for iterating through minimum coverage statistics
    unsigned int max_min_coverage = 0; // will be adjusted below.
    unsigned int cov_step = 1; 
    
    float max_con_ratio = 0.0;
    float con_step = 0.05;
    
    float max_min_entropy = 0.0;
    float entropy_step = 0.1;
    
    bool WRITE_COVERAGE = false;
    string COVERAGE_OUTPUT_FILENAME;
    
    
    // error-containing kmer removal:
    bool prune_error_kmers = true;
    float min_ratio_non_error = 0.005f;
    

    bool PARALLEL_IWORM = false;
    bool SINGLE_PHASE = false;

    
    bool SUPER_READS_MODE = false;
    
    /***************************************/
    /***   Option Processing   ************/
    /**************************************/
    
        
    try {
        
        ArgProcessor args(argc, argv);
        
        /* Check for essential options */
        if (args.isArgSet("--help") || args.isArgSet("-h") 
            ||  (! (args.isArgSet("--reads") || args.isArgSet("--kmers") || args.isArgSet("--kmer_files_listing") ) )
            
            || (! (args.isArgSet("--run_inchworm") 
                   || args.isArgSet("--make_super_reads")
                   || args.isArgSet("--checkFastaPath")
                   || args.isArgSet("--threadFasta")
                   || args.isArgSet("--describe_kmers")
                   
                   ) 
                )
            ) {
            cerr  << usage(args) << endl << endl;
            
            return 1;
        }
        
        /* Process option values */

        // ************************** //
        //  NOTE:  commenting out some of the options that were experimental and now just clutter the code // bhaas, 1-2012
        // ************************* //

        
        if (args.isArgSet("--reads"))
        		reads_fasta_filename = args.getStringVal("--reads");
        
        if (args.isArgSet("--kmers"))
	        kmers_fasta_filename = args.getStringVal("--kmers");

        if (args.isArgSet("--kmer_files_listing"))
            kmer_files_list_filename = args.getStringVal("--kmer_files_listing");
        

        if (args.isArgSet("-K")) {
            kmer_length = args.getIntVal("-K"); 
            cerr << "Kmer length set to: " << kmer_length << endl;
        }
        
        if (args.isArgSet("--minKmerCount")) {
            
            min_kmer_count = args.getIntVal("--minKmerCount");
            cerr << "Min coverage set to: " << min_kmer_count << endl;
            
        }
        
        if (args.isArgSet("-L")) {
            MIN_ASSEMBLY_LENGTH = args.getIntVal("-L");
            cerr << "Min assembly length set to: " << MIN_ASSEMBLY_LENGTH << endl; 
        }
        if (args.isArgSet("--min_assembly_coverage")) {
            MIN_ASSEMBLY_COVERAGE = args.getIntVal("--min_assembly_coverage");
            cerr << "Min assembly coverage set to: " << MIN_ASSEMBLY_COVERAGE << endl;
        }

        if (args.isArgSet("--make_super_reads")) {
            SUPER_READS_MODE = true;
            run_inchworm = true;
            cerr << "Running Inchworm in Super-Reads mode" << endl;
        }
        
        /* CRAWL was experimental
           
        if (args.isArgSet("--crawl")) {
            CRAWL = true;
            cerr << "Crawling enabled" << endl;
            if (args.isArgSet("--crawl_length")) {
                crawl_length = args.getIntVal("--crawl_length");
                cerr <<  "Crawl length set to: " << crawl_length << endl;
            }
        }
        */
        
        /*   pacman was experimental
             
        if (args.isArgSet("--pacman")) {
            PACMAN = true;
            cerr << "PACMAN enabled." << endl;
        }
        
        */

        if (args.isArgSet("--monitor")) {
            IRKE_COMMON::MONITOR = args.getIntVal("--monitor");
            cerr << "Monitor turned on, set to: " << IRKE_COMMON::MONITOR << endl;
        }

        if (args.isArgSet("--keep_tmp_files")) {
            IRKE_COMMON::KEEP_TMP_FILES = true;
            cerr << "-retaining tmp files" << endl;
        }
        
        if (args.isArgSet("--min_con_ratio")) {
            MIN_CONNECTIVITY_RATIO = args.getFloatVal("--min_con_ratio");
        }
        
        if (args.isArgSet("--DS")) {
            double_stranded_flag = true;
            cerr << "double stranded mode set" << endl;
        }
        
        if (args.isArgSet("--min_seed_entropy")) {
            MIN_SEED_ENTROPY = args.getFloatVal("--min_seed_entropy");
            cerr << "Min seed entropy set to: " << MIN_SEED_ENTROPY << endl;
        }
        if (args.isArgSet("--min_seed_coverage")) {
            MIN_SEED_COVERAGE = args.getIntVal("--min_seed_coverage");
            cerr << "min seed coverage set to: " << MIN_SEED_COVERAGE << endl;
        }
        
        if (args.isArgSet("--min_any_entropy")) {
            min_any_entropy = args.getFloatVal("--min_any_entropy");
            cerr << "min entropy set to: " << min_any_entropy << endl;
        }
        
        if (args.isArgSet("--run_inchworm")) {
            run_inchworm = true;
        
            /*
            if (args.isArgSet("-R")) {
                MAX_RECURSION = args.getIntVal("-R");
                cerr << "Setting max recursion to: " << MAX_RECURSION << endl;
            }
            */
            
        }


        // multithreading settings:
        omp_set_dynamic(false);
        if (args.isArgSet("--num_threads")) {
            IRKE_COMMON::NUM_THREADS = args.getIntVal("--num_threads");
            omp_set_num_threads(IRKE_COMMON::NUM_THREADS);
            cerr << "setting number of threads to: " << IRKE_COMMON::NUM_THREADS << endl;
        }
        else {
            IRKE_COMMON::NUM_THREADS = omp_get_max_threads();
        }
        if (IRKE_COMMON::MAX_THREADS_READ_PARSING > IRKE_COMMON::NUM_THREADS) {
            IRKE_COMMON::MAX_THREADS_READ_PARSING = IRKE_COMMON::NUM_THREADS;
        }
        
        
        if (args.isArgSet("--checkFastaPath")) {
            cds_fasta_filename = args.getStringVal("--checkFastaPath");
            search_CDS_flag = true;
            
            if (args.isArgSet("--max_min_cov")) {
                max_min_coverage = args.getIntVal("--max_min_cov");
                
                if (args.isArgSet("--cov_step")) {
                    cov_step = args.getIntVal("--cov_step");
                }
            }
            
            /*
            if (args.isArgSet("--max_con_ratio")) {
                max_con_ratio = args.getFloatVal("--max_con_ratio");
                
                if (args.isArgSet("--con_step")) {
                    con_step = args.getFloatVal("--con_step");
                }
            }
            */
            
            if (args.isArgSet("--max_min_entropy")) {
                max_min_entropy = args.getFloatVal("--max_min_entropy");
                
                if (args.isArgSet("--entropy_step")) {
                    entropy_step = args.getFloatVal("--entropy_step");
                }
                
            }
        }
        
        
        if (args.isArgSet("--threadFasta")) {
            cds_fasta_filename = args.getStringVal("--threadFasta");
            thread_CDS_flag = true;
        }
        
       
        if (args.isArgSet("--describe_kmers")) {
            describe_kmers = true;
        }
       
        
        if (args.isArgSet("--coverage_outfile")) {
            WRITE_COVERAGE = true;
            COVERAGE_OUTPUT_FILENAME = args.getStringVal("--coverage_outfile");
        }
        
        /*
        if (args.isArgSet("--reassembleIworm")) {
            reassembleIworm = true;
            cerr << "Running in re-assembly mode." << endl;
        }
        
       

        if (args.isArgSet("--pruneSingletonKmersReadInterval")) {
            pruneSingletonKmersReadInterval = args.getLongVal("--pruneSingletonKmersReadInterval");
        }
        
        */
        
 
        // kmer error removal options
        if (args.isArgSet("--no_prune_error_kmers")) {
            prune_error_kmers = false;
        }
        
        if (prune_error_kmers && args.isArgSet("--min_ratio_non_error")) {
            min_ratio_non_error = args.getFloatVal("--min_ratio_non_error");
            cerr << "Set to prune kmers below min ratio non-erro: " << min_ratio_non_error << endl;
        }

        if (args.isArgSet("--PARALLEL_IWORM")) {
            PARALLEL_IWORM = true;
            IRKE_COMMON::__DEVEL_no_kmer_sort = true;  // not sorting the kmers. //FIXME: this needs to be less hacky
            cerr << "-setting parallel iworm mode." << endl;
            
            if (args.isArgSet("--SINGLE_PHASE")) {
                SINGLE_PHASE = true;
                cerr << "-setting single phase parallel iworm build." << endl;
            }
        }
                
        // process developer parameters
        if (args.isArgSet("--DEVEL_no_kmer_sort")) {
            IRKE_COMMON::__DEVEL_no_kmer_sort = true;
            cerr << "Setting IRKE_COMMON::__DEVEL_no_kmer_sort = true" << endl;
        }
        if (args.isArgSet("--DEVEL_no_greedy_extend")) {
            IRKE_COMMON::__DEVEL_no_greedy_extend = true;
            cerr << "Setting IRKE_COMMON::__DEVEL_no_greedy_extend = true" << endl;
        }
        if (args.isArgSet("--DEVEL_no_tie_breaking")) {
            IRKE_COMMON::__DEVEL_no_tie_breaking = true;
            cerr << "Setting IRKE_COMMON::__DEVEL_no_tie_breaking = true" << endl;
        }
        if (args.isArgSet("--DEVEL_zero_kmer_on_use")) {
            IRKE_COMMON::__DEVEL_zero_kmer_on_use = true;
            cerr << "Setting IRKE_COMMON::__DEVEL_zero_kmer_on_use = true" << endl;
        }
        
        if (args.isArgSet("--DEVEL_rand_fracture_prob")) {
            IRKE_COMMON::__DEVEL_rand_fracture_prob = args.getFloatVal("--DEVEL_rand_fracture_prob");
            IRKE_COMMON::__DEVEL_rand_fracture = true;
            cerr << "Setting --DEVEL_rand_fracture_prob to " << IRKE_COMMON::__DEVEL_rand_fracture_prob << endl;
        }
        
        
        /* adjust values of max if not set directly or set improperly */
        if (max_min_coverage < min_kmer_count) {
            max_min_coverage = min_kmer_count;
        }
        
        if (max_con_ratio < MIN_CONNECTIVITY_RATIO) {
            max_con_ratio = MIN_CONNECTIVITY_RATIO;
        }
        
        if (max_min_entropy < min_any_entropy) {
            max_min_entropy = min_any_entropy;
        }
                
    }
    
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    
    /**********************************************/
    /********  Run Inchworm  **********************/
    /**********************************************/
    
    
    IRKE irke(kmer_length, MAX_RECURSION, MIN_SEED_ENTROPY, MIN_SEED_COVERAGE, min_any_entropy, PACMAN, CRAWL, crawl_length, double_stranded_flag);
    // note PACMAN & CRAWL were experimental, didn't seem to help, and should be stripped out at some point.  (bhaas)


    if (SINGLE_PHASE) {
        irke.TWO_PHASE = false;
    }
    
    if (pruneSingletonKmersReadInterval > 0) {
        irke.set_prune_singleton_read_interval(pruneSingletonKmersReadInterval);
    }
    
    
    // Kmer catalog construction
    if (reads_fasta_filename.length()) {
        // let inchworm do it (pretty slow)
    		irke.build_graph(reads_fasta_filename, reassembleIworm, false);

    } 
    else if (kmers_fasta_filename.length()) {
    	// read an existeing kmer catalog (such as from Meryl)
        irke.build_graph(kmers_fasta_filename, reassembleIworm, true);  
        // The  reassembleIworm idea (from different k-mer lengths) didn't pan out. Should remove it. (bhaas)
    }
    else if (kmer_files_list_filename.length()) {
        vector<string> kmer_files_list = parse_file_listing(kmer_files_list_filename);
        for (unsigned int i = 0; i < kmer_files_list.size(); i++) {
            string kmer_filename = kmer_files_list[i];
            cerr << endl << endl << "// parsing kmer file: " << kmer_filename << endl;
            irke.build_graph(kmer_filename, reassembleIworm, true);
        }
    }
    
    
    if (search_CDS_flag) { // define the 'oracle' transcript set
        
        // iterate over min coverage values
        for (unsigned int minC = min_kmer_count; minC <= max_min_coverage; minC += cov_step) {
            
            // iterate over min connectivity values
            for (float min_conn = MIN_CONNECTIVITY_RATIO; min_conn <= max_con_ratio; min_conn += con_step) {
                
                for (float min_entropy = min_any_entropy; min_entropy <= max_min_entropy; min_entropy += entropy_step) {
                    
                    examine_CDS_paths(irke, cds_fasta_filename, minC, min_conn, min_entropy, WRITE_COVERAGE, COVERAGE_OUTPUT_FILENAME);
                    
                }
                
            }
        }
    }
    
    else if (thread_CDS_flag) {
        
        thread_sequences_through_graph(irke, cds_fasta_filename);

        
    }
    else {
        
                
        if (min_kmer_count > 1 || min_any_entropy > 0 || prune_error_kmers) {
            cerr << "Pruning kmers (min_kmer_count=" << min_kmer_count
                 << " min_any_entropy=" << min_any_entropy
                 << " min_ratio_non_error=" << min_ratio_non_error << ")"<< endl;
            int start = time(NULL);
            irke.prune_some_kmers(min_kmer_count, min_any_entropy, prune_error_kmers, min_ratio_non_error);
            int end = time(NULL);
            int pruning_time = end-start;
            cerr << "\tPruning time: " << pruning_time << " seconds = " << (float)pruning_time/60.0 << " minutes." << endl;  
        
            cerr << endl << "TIMING PRUNING " << pruning_time << " s." << endl;
        }

        
        if (SUPER_READS_MODE) {
            cerr << "Pruning branched kmers under super-reads mode:" << endl;
            irke.prune_branched_kmers();
        }
        
        if (describe_kmers) {
            irke.describe_kmers();
        }
    
            
        if (run_inchworm) { // note this would have to be true since either or both are required.
            
            // sort kmers by abundance descendingly (unless PARALLEL_IWORM option, in which case the order is just the hash iterator order)
        
            cerr << "-populating the kmer seed candidate list." << endl;
            irke.populate_sorted_kmers_list(); 

            cerr << "-beginning inchworm contig assembly." << endl;
            int start = time(NULL);
            irke.compute_sequence_assemblies(MIN_CONNECTIVITY_RATIO, MIN_ASSEMBLY_LENGTH, MIN_ASSEMBLY_COVERAGE, 
                                             PARALLEL_IWORM, WRITE_COVERAGE, COVERAGE_OUTPUT_FILENAME);
            
            int end = time(NULL);
            int iworm_assembly_time = end-start;
            cerr << "\tIworm contig assembly time: " << iworm_assembly_time << " seconds = " 
                 << (float)iworm_assembly_time/60.0 << " minutes." << endl;
            
            cerr << endl << "TIMING CONTIG_BUILDING " << iworm_assembly_time << " s." << endl;
        }
        
    }
    

    int prog_end = time(NULL);

    int prog_runtime = prog_end - prog_start;

    cerr << endl << "TIMING PROG_RUNTIME " << prog_runtime << " s." << endl;
    
    
    return (0);
    
}



void examine_CDS_paths(IRKE& irke, string cds_fasta_filename, unsigned int min_cov, float min_connectivity, float min_entropy,
					   bool WRITE_COVERAGE, string COVERAGE_OUTPUT_FILENAME) {
    
    
    Fasta_reader fasta_reader(cds_fasta_filename);
    
    ofstream coverage_writer;
    if (WRITE_COVERAGE) {
        coverage_writer.open(COVERAGE_OUTPUT_FILENAME.c_str());
        if (! coverage_writer.is_open()) {
            throw(stacktrace() + "Error, cannot write to file: " + COVERAGE_OUTPUT_FILENAME);
        }
    }
    
    
    while (fasta_reader.hasNext()) {
        
        Fasta_entry fe = fasta_reader.getNext();
        
        string accession = fe.get_accession();
        
        string seq = fe.get_sequence();
        
        vector<unsigned int> coverage_counter;
        
        cout << accession << "\tCov: " << min_cov << "\tCon: " << min_connectivity << "\tE: " << min_entropy;
        
        if (irke.sequence_path_exists(seq, min_cov, min_entropy, min_connectivity, coverage_counter)) {
            
            cout << "\tT" << endl;
            
        }
        else {
            
            cout << "\tF" << endl;
        }
        
        if (WRITE_COVERAGE) {
            
            coverage_writer << ">" << accession << endl;
            
            for (unsigned int i = 0; i < coverage_counter.size(); i++) {
                coverage_writer << coverage_counter[i];
                if ( (i +1) % 30 == 0) {
                    coverage_writer << endl;
                }
                else {
                    coverage_writer << " ";
                }
            }
            coverage_writer << endl;
        }
        
        
    }
    
    if (WRITE_COVERAGE) {
        coverage_writer.close();
    }
    
    
    return;
}



void thread_sequences_through_graph(IRKE& irke, string cds_fasta_filename) {
    
    Fasta_reader fasta_reader(cds_fasta_filename);
    
    while (fasta_reader.hasNext()) {
        
        Fasta_entry fe = fasta_reader.getNext();
        
        string accession = fe.get_accession();
        
        string seq = fe.get_sequence();
        
        cout << "// " << accession << endl;
        
        cout << irke.thread_sequence_through_graph(seq) << endl << endl;
        
    }
    
    return;
}

string usage (ArgProcessor args) {
    
    bool show_advanced = args.isArgSet("--show_advanced");
    
    stringstream usage_info;
    
    usage_info 
        << endl << endl
        << "Usage:" << endl
        << "  inchworm --reads <filename.fasta> --run_inchworm [opts] " << endl
        << "  inchworm --kmers <filename.fasta> --run_inchworm [opts] " << endl 
        << "  inchworm --kmer_files_listing <kmer_file_list.txt> --run_inchworm [opts] " << endl << endl
        

        << "**Required" << endl
        << "  --reads  <str>             " << ":fasta file containing reads" << endl
        << "  --kmers  <str>             " << ":fasta file containing kmers" << endl
        << "  --kmer_files_listing <str> " << ":file listing filenames containing kmers" << endl
        << endl;
    
    
    
    usage_info
        << endl
        << "** Run modes:" << endl
        << "  --run_inchworm           " << ":run inchworm, report sequences" << endl
        << "  --make_super_reads       " << ":generates super-reads" << endl;
    
    if (show_advanced) {
        usage_info 
            << "    * advanced: " << endl
            << "    --describe_kmers         " << ":just describe kmers in the graph" << endl
            << "    --checkFastaPath <str>   " << ":check fastea seqs for path in kmer graph" << endl
            << "    --threadFasta <str>      " << ":thread sequences through the kmer graph, report kmer" << endl;
        
    }
    
    
    usage_info
        << endl
        << "** General opts" << endl
        << "  -K <int>                         " << ":kmer length (default: 25, meaning 24mer overlaps)  max = 32 (stored as 64-bit integers, 2-base char encoding)" << endl
        << "  --minKmerCount <int>             " << ":min kmer count, lower abundant kmers are pruned from graph (default: 1)" << endl
        << "  -L <int>                         " << ":min contig length to be reported (default: 25)" << endl
        << "  --min_assembly_coverage  <int>   " << ":min kmer coverage of an assembled contig to be reported (default: 2)" << endl
        << "  --coverage_outfile <str>         " << ":file to store per-base kmer coverage for contigs" << endl
        << "  --DS                             " << ":double-stranded RNA-Seq mode (not strand-specific)" << endl
        << "  --no_prune_error_kmers           " << ":disable pruning of kmers that occur at below --min_ratio_non_error " << endl
        << "  --min_ratio_non_error <float>    " << ":min ratio for low/high alternative extension that is not an error (default: 0.005)" << endl
        << "  --num_threads <int>              " << ":number of threads to use. (by default, uses OMP_NUM_THREADS env var setting)" << endl
        << "  --PARALLEL_IWORM                 " << ":run the contig building in parallel; by default, only does parallel read parsing." << endl
      
        ;
    
    if (show_advanced) {
        usage_info 
            << "    * advanced opts:" << endl
            // << "    -R <int>                  " << ":maximum recursion length (default: 1)" << endl
            << "    --min_any_entropy <float> " << ":kmer must have at least this entropy to be added to the graph (default: 0)" << endl
            << "    --min_seed_entropy <float> " << ":min seed kmer entropy (default: 1.5)" << endl
            << "    --min_seed_coverage <int> " << ":min seed kmer coverage (default: 2)" << endl
            
            // << "    --crawl                   " << ":drag tail of inchworm by --crawl_length" << endl
            // << "    --crawl_length  <int>     " << ":length of tail end crawl (default: 1)" << endl
            // << "    --pacman                  " << ":exclude kmers as they're visited in descending order of counts" << endl
            // << "    --pruneSingletonKmersReadInterval <long> " << ":at num reads parsed interval, prune singleton kmers (default: 0)" << endl
            // << "    --min_con_ratio <float> " << ":minimum connectivity ratio (min/max) between neighboring kmers (default: 0)" << endl
            // << "    --reassembleIworm       " << ":redo assembly using previous inchworm outputs. Post-process to vary-K" << endl;

            << "    * developer flags: " << endl
            << "    --DEVEL_no_kmer_sort      "  << ":no sorting of kmer by abundance values occurs. Kmer seeds chosen according to hash iterator order." << endl
            << "    --DEVEL_no_greedy_extend  "  << ":instead of greedy contig extensions, any non-zero kmer is chosen to extend." << endl
            << "    --DEVEL_no_tie_breaking   "  << ":if tie, choose either path randomly" << endl
            << "    --DEVEL_zero_kmer_on_use  "  << ":zero kmer counts on use rather than after contig construction" << endl
            << "    --DEVEL_rand_fracture_prob <float>  " << " :probability of fracturing an inchworm contig at each extension (default: 0 (off)" << endl
            
            << "    --SINGLE_PHASE               " << ":in parallel mode, do a single iworm construction phase instead of default 2-phase" << endl
            << "    --keep_tmp_files             " << ":retain temporary files." << endl    
            

            ; 

        
    }
    
    
    usage_info 
        << endl
        << "** Misc opts." << endl
        << "  --help                   " << ":produce this menu." << endl
        << "  --monitor <int>          " << ":verbosity. ( '1': recommended, '2': for debugging ) " << endl
        << "  --show_advanced          " << ":more advanced options (mostly experimental)" << endl 
        
        << endl << endl << endl;
    
    
    return(usage_info.str());
	
}


 vector<string> parse_file_listing (string filename) {

     ifstream in;
     in.open(filename.c_str());
     if (! in.is_open()) {
         stringstream errstr;
         errstr << "Error, cannot open file: " << filename;
         throw(errstr.str());
     }

     string line;
     vector<string> lines;
     while (! in.eof()) {
         if (getline(in, line)) {
             lines.push_back(line);
         }
     }
     in.close();

     return(lines);
 }

         
     
     
