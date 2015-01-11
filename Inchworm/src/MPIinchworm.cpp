#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <unistd.h> 

#include "irke_common.hpp"
#include "argProcessor.hpp"
#include "Fasta_entry.hpp"
#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"

#include <mpi.h>
#include <omp.h>


using namespace std;

// confiuration of globals
bool IRKE_COMMON::__DEVEL_no_kmer_sort = true;  // do not sort kmers in parallel iworm mode.
int IRKE_COMMON::NUM_THREADS = 1;
int IRKE_COMMON::MAX_THREADS_READ_PARSING = 6;
unsigned int IRKE_COMMON::MONITOR = 0;

unsigned int MAX_TEST_KMERS = 10;

unsigned int NUM_MPI_NODES;
unsigned int THIS_MPI_NODE_ID;
bool MASTER_SLAVE_MODE = false;
bool WRITE_KMER_FILES = false;
bool MONITOR_MPI_COMMUNICATION = false;
bool KEEP_TMP_FILES = false;
int NUM_MPI_CONTIG_BUILDERS = -1; // turn off by default. All nodes are both contig builders and kmer servers.

unsigned int KMER_SERVER_PHASE = 0;


float MIN_CONNECTIVITY_RATIO = 0.0;
unsigned int MIN_ASSEMBLY_LENGTH; // minimum length of an inchworm assembly for reporting.
unsigned int MIN_ASSEMBLY_COVERAGE = 2; // minimum average kmer coverage for assembly to be reported.
float MIN_SEED_ENTROPY = 1.5; // minimum entropy for a Kmer to seed in inchworm assembly construction.
unsigned int MIN_SEED_COVERAGE = 2; // minimum kmer coverage for a seed.
bool DOUBLE_STRANDED_MODE = false;  // strand-specific by default

bool THIS_NODE_DONE = false;
bool ALL_NODES_DONE = false;

bool MPI_SLEEP_INTERVAL = 5; // seconds between checks for completion status at end.

string TOKEN;


enum { // work requests
    GENERAL_WORK_REQUEST, GET_RIGHT_GREEDY_EXTENSION, GET_LEFT_GREEDY_EXTENSION, SET_ALL_NODES_DONE, GET_FINISHED_STATUS, CLEAR_KMER, ADD_KMER,
    
    // responses
    GENERAL_WORK_RESPONSE, SENDING_RIGHT_GREEDY_EXTENSION, SENDING_LEFT_GREEDY_EXTENSION, SENDING_FINISHED_STATUS, SENDING_KMER_ADD_CONFIRM,

    // other
    ERROR_ENCOUNTERED
};

map<unsigned int,string> enum_description;



// function prototypes
string usage(ArgProcessor args);
kmer_int_type_t get_central_kmer(kmer_int_type_t kmer, unsigned int kmer_length);
kmer_int_type_t get_central_right_kmer(kmer_int_type_t kmer, unsigned int kmer_length);
kmer_int_type_t get_central_left_kmer(kmer_int_type_t kmer, unsigned int kmer_length);

unsigned int get_node_for_central_kmer (kmer_int_type_t central_kmer, unsigned int kmer_length);

vector<kmer_int_type_t>  build_inchworm_contig_from_seed(kmer_int_type_t kmer, KmerCounter& kcounter, float min_connectivity, unsigned int& total_counts);
vector<kmer_int_type_t> join_forward_n_reverse_paths(vector<kmer_int_type_t>& reverse_path, kmer_int_type_t seed_kmer_val, vector<kmer_int_type_t>& forward_path);
vector<kmer_int_type_t> inchworm (KmerCounter& kcounter, char direction, kmer_int_type_t kmer, Kmer_visitor& visitor, float min_connectivity);

bool is_good_seed_kmer(KmerCounter& kcounter, kmer_int_type_t kmer, unsigned int kmer_count, unsigned int kmer_length, float min_connectivity);
string reconstruct_path_sequence(KmerCounter& kcounter, vector<kmer_int_type_t>& path);

void zap_kmers(KmerCounter& kcounter, vector<kmer_int_type_t>& kmer_path);
void run_MPI_kmer_server(KmerCounter& kcounter);
void run_MPI_master_all_completion_check(int phase_val);

void test_MPI(KmerCounter& kcounter);
kmer_int_type_t get_greedy_extension_from_MPI_node(unsigned int node_for_central_kmer, kmer_int_type_t kmer, char direction, unsigned int kmer_length);

string get_MPI_node_filename(unsigned int node_id);

string add_fasta_seq_line_breaks(string& sequence, int interval);
kmer_int_type_t extract_best_seed(vector<kmer_int_type_t>& kmer_vec, KmerCounter& kcounter, float min_connectivity);
void add_kmer_to_kcounter(KmerCounter& kcounter, kmer_int_type_t& kmer, unsigned int& count, ofstream& logstream);


int main (int argc, char* argv[]) {

  stringstream timings_text;

  time_t prog_start_time;
  time(&prog_start_time);
  

    // ----------------------------
    // Begin MPI configuration: (adapted from some example code from somewhere...)
    //-----------------------------
    
    int tid,nthreads;
    char *cpu_name;
    double time_initial,time_current,time_used;

    

    /* add in MPI startup routines */
    /* 1st: launch the MPI processes on each node */
    int provided = -1;
    MPI_Init_thread(&argc,&argv, MPI_THREAD_MULTIPLE, &provided);

    /* 2nd: request a thread id, sometimes called a "rank" from
          the MPI master process, which has rank or tid == 0
    */
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);
    THIS_MPI_NODE_ID = tid; // set my global var.
    

    if (provided != MPI_THREAD_MULTIPLE) {
      if (THIS_MPI_NODE_ID == 0) {
	cerr << "MPIinchworm cannot run without MPI_THREAD_MULTIPLE support enabled." << endl
	     << "\tmax support set to: " << provided << ", but need: " << MPI_THREAD_MULTIPLE << endl;
      }
      exit(1);
    }

    time_initial  = MPI_Wtime();
    
    /* 3rd: this is often useful, get the number of threads
          or processes launched by MPI (max(tid) +1)
    */
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
    NUM_MPI_NODES = nthreads; // set my global var

    // needed for nested parallelism in contig building.
    omp_set_nested(1);
    
    // populate the enum descriptor
    enum_description[GENERAL_WORK_REQUEST] = "GENERAL_WORK_REQUEST";
    enum_description[GET_RIGHT_GREEDY_EXTENSION] = "GET_RIGHT_GREEDY_EXTENSION";
    enum_description[GET_LEFT_GREEDY_EXTENSION] = "GET_LEFT_GREEDY_EXTENSION";
    enum_description[SET_ALL_NODES_DONE] = "SET_ALL_NODES_DONE";
    enum_description[GET_FINISHED_STATUS] = "GET_FINISHED_STATUS";
    enum_description[ADD_KMER] = "ADD_KMER";
    enum_description[CLEAR_KMER] = "CLEAR_KMER";
    enum_description[GENERAL_WORK_RESPONSE] = "GENERAL_WORK_RESPONSE";
    enum_description[SENDING_RIGHT_GREEDY_EXTENSION] = "SENDING_RIGHT_GREEDY_EXTENSION";
    enum_description[SENDING_LEFT_GREEDY_EXTENSION] = "SENDING_LEFT_GREEDY_EXTENSION";
    enum_description[SENDING_FINISHED_STATUS] = "SENDING_FINISHED_STATUS";
    enum_description[SENDING_KMER_ADD_CONFIRM] = "SENDING_KMER_ADD_CONFIRM";
    enum_description[ERROR_ENCOUNTERED] = "ERROR_ENCOUNTERED";
    
    string fasta_filename; // single file of k-mers
    unsigned int kmer_length = 25;
    unsigned int min_kmer_count = 1;
    MIN_ASSEMBLY_LENGTH = kmer_length;
    
    int READ_TYPE;
    enum { KMER, READ };


    // error-containing kmer removal:
    bool prune_error_kmers = true;
    float min_ratio_non_error = 0.05f;
    
    
    // argument processing.
    ArgProcessor args(argc, argv);
    
    /* Check for essential options */
    if (args.isArgSet("--help") || args.isArgSet("-h") 
        ||  (! (args.isArgSet("--kmers") || args.isArgSet("--reads") ) )
                 || (! args.isArgSet("--token") )
        ) {
        
        cerr  << usage(args) << endl << endl;
        
        return 1;
    }
        
    
    // required params
    if (args.isArgSet("--kmers")) {
        fasta_filename = args.getStringVal("--kmers");
        READ_TYPE = KMER;
    }
    else if (args.isArgSet("--reads")) {
        fasta_filename = args.getStringVal("--reads");
        READ_TYPE = READ;
    }
    else {
        cerr << "Error, must specify --kmers or --reads";
        exit(4);
    }
    
    TOKEN = args.getStringVal("--token");
    
    // optional args
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
    
    if (args.isArgSet("--monitor")) {
        IRKE_COMMON::MONITOR = args.getIntVal("--monitor");
        cerr << "Monitor turned on, set to: " << IRKE_COMMON::MONITOR << endl;
    }
        
    if (args.isArgSet("--min_con_ratio")) {
        MIN_CONNECTIVITY_RATIO = args.getFloatVal("--min_con_ratio");
    }
        
    if (args.isArgSet("--DS")) {
        DOUBLE_STRANDED_MODE = true;
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

    
    // some testing parameters.
    if (args.isArgSet("--max_test_kmers")) {
        MAX_TEST_KMERS = args.getIntVal("--max_test_kmers");
        cerr << "-SETTING MAX_TEST_KMERS to " << MAX_TEST_KMERS << endl;
    }
    if (args.isArgSet("--master_slave_mode")) {
        MASTER_SLAVE_MODE = true;
        cerr << "-MASTER_SLAVE_MODE in effect. " << endl;
    }
    if (args.isArgSet("--num_mpi_contig_builders")) {
       
        if (MASTER_SLAVE_MODE) {
            cerr << "ERROR, cannot be in master/slave mode and setting --num_mpi_contig_builders, conflicting run modes. " << endl;
            exit(1);
        }
        
        NUM_MPI_CONTIG_BUILDERS = args.getIntVal("--num_mpi_contig_builders");
        if (NUM_MPI_CONTIG_BUILDERS < 1) {
            cerr << "Error, num_mpi_contig_builders to be set to " << NUM_MPI_CONTIG_BUILDERS << ", not allowed, must be a positive integer. ";
            exit(1);
        }
        
        if (THIS_MPI_NODE_ID == 0) 
            cerr << "--setting NUM_MPI_CONTIG_BUILDERS to " << NUM_MPI_CONTIG_BUILDERS << endl;
        
    }
    

    if (args.isArgSet("--write_kmer_files")) {
        cerr << "-node will write its kmer database file" << endl;
        WRITE_KMER_FILES = true;
    }
    if (args.isArgSet("--monitor_MPI")) {
        cerr << "-node is monitoring MPI communications." << endl;
        MONITOR_MPI_COMMUNICATION = true;
    }
    if (args.isArgSet("--keep_tmp_files")) {
        KEEP_TMP_FILES = true;
    }
    // end of testing params
    
    
    if (args.isArgSet("--num_threads")) {
        IRKE_COMMON::NUM_THREADS = args.getIntVal("--num_threads");
    }
    omp_set_num_threads(IRKE_COMMON::NUM_THREADS);
    cerr << "setting number of threads to: " << IRKE_COMMON::NUM_THREADS << endl;
    
    if (IRKE_COMMON::MAX_THREADS_READ_PARSING > IRKE_COMMON::NUM_THREADS) {
        IRKE_COMMON::MAX_THREADS_READ_PARSING = IRKE_COMMON::NUM_THREADS;
    }
    
    
    // kmer error removal options
    if (args.isArgSet("--no_prune_error_kmers")) {
        prune_error_kmers = false;
    }
    if (prune_error_kmers && args.isArgSet("--min_ratio_non_error")) {
        min_ratio_non_error = args.getFloatVal("--min_ratio_non_error");
        cerr << "Set to prune kmers below min ratio non-erro: " << min_ratio_non_error << endl;
    }
    

    

    if (NUM_MPI_CONTIG_BUILDERS > 0 && NUM_MPI_CONTIG_BUILDERS == NUM_MPI_NODES) {
        cerr << "WARNING:  All nodes are both kmer servers and contig builders. This is the default setting, so no need to parameterize it as such. " << endl;
        
    }

    
    cerr << "Number of MPI nodes: " << nthreads << endl;
    
    /* on EVERY process, allocate space for the machine name */
    cpu_name    = (char *)calloc(80,sizeof(char));

    /* get the machine name of this particular host ... well
       at least the first 80 characters of it ... */
    gethostname(cpu_name,80);
    time_current  = MPI_Wtime();
    time_used  = time_current - time_initial;
    fprintf(stderr, "%.3f tid=%i : hello MPI user: machine=%s [NCPU=%i]\n",
           time_used, tid, cpu_name, nthreads);


    //--------------------------------------------------------------------
    // Determine if this node is to be a builder or kmer server (or both)
    //--------------------------------------------------------------------

    bool is_assembly_node = true; 
        
    if (MASTER_SLAVE_MODE) {
        if (THIS_MPI_NODE_ID == 0) {
            is_assembly_node = true;
        }
        else {
            is_assembly_node = false;
        }
    }
    else if (NUM_MPI_CONTIG_BUILDERS > 0) {

        if (THIS_MPI_NODE_ID < NUM_MPI_CONTIG_BUILDERS) {
            is_assembly_node = true;
        }
        else {
            is_assembly_node = false;
        }
        
    }
    
    //-----------------------------------------------------
    // main data store for node
    KmerCounter kcounter(kmer_length, DOUBLE_STRANDED_MODE); 
    //------------------------------------------------------

    
    time_t kmer_db_build_start_time;
    time(&kmer_db_build_start_time);
    
    #pragma omp parallel sections num_threads(2)
    {
        
        #pragma omp section
        {             
            // everyone opens a communication thread, regardless of whether it's storing kmers or not.
            
	  if (NUM_MPI_NODES == 1) {
	    // no reason to run kmer server
	    cerr << "** Phase 1: Only 1 MPI node, no reason to do MPI communication. exiting run_MPI_kmer_server()" << endl;
	   
	  }
	  else {
	    run_MPI_kmer_server(kcounter);
            
	  }
	}
        
        #pragma omp section 
        {
            // everyone participates in reading a part of the kmer file
            // and getting the kmer servers populated.
            
            // figure out which section of the kmer fasta file we're supposed to use:
            long file_length = 0; // init, set below
            long this_mpi_section_start = 0;
            long this_mpi_section_end = -1;
            
            if (NUM_MPI_NODES > 1) {
                
                
                ifstream fasta_file_reader (fasta_filename.c_str());
                fasta_file_reader.seekg(0, fasta_file_reader.end);
                file_length = fasta_file_reader.tellg();
                fasta_file_reader.seekg(0, fasta_file_reader.end);
                fasta_file_reader.close();
                
                
                long mpi_section_length = file_length / NUM_MPI_NODES;
                this_mpi_section_start = THIS_MPI_NODE_ID * mpi_section_length;
                this_mpi_section_end = this_mpi_section_start + mpi_section_length;
            }
            
            //---------------------------------------------
            // Kmer partitioning among the nodes: 
            // Populate kmer hashtable on each node
            // where each node gets a subset of the kmers
            //---------------------------------------------
            
            
            Fasta_reader fasta_reader(fasta_filename, this_mpi_section_start, this_mpi_section_end);
            
            stringstream filenamebuilder;
            filenamebuilder << "tmp." << TOKEN << ".kmers.tid_" << tid;
            ofstream filewriter;
            
            
            if (WRITE_KMER_FILES) {
                filewriter.open(filenamebuilder.str().c_str());
            }
            
            
            filenamebuilder.str("");
            filenamebuilder << "mpi_kmer_file_reader.mpi_" << THIS_MPI_NODE_ID << ".log";
            ofstream kmerReaderLog;
            if (MONITOR_MPI_COMMUNICATION)
                kmerReaderLog.open(filenamebuilder.str().c_str());
            
            MPI_Status MPI_status;
            
            unsigned long kmer_counter = 0;
            while (true) {

                if (! fasta_reader.hasNext())
                    break;

                Fasta_entry fe = fasta_reader.getNext();
                string seq = fe.get_sequence();
                
                if (seq == "") continue;
                
                if (seq.length() < kmer_length) {
                    continue;
                }
                
                kmer_counter++;

				unsigned int count = 1;
                if (READ_TYPE == KMER)
                    count = atoi(fe.get_header().c_str());
                
                for (int i = 0; i <= seq.length() - kmer_length; i++) {
                    
                    string kmer_s = seq.substr(i, kmer_length);

                    if (contains_non_gatc(kmer_s))
                        continue;

                    kmer_int_type_t kmer = kcounter.get_kmer_intval(kmer_s);
                    kmer_int_type_t central_kmer = get_central_kmer(kmer, kmer_length);
                
                
                    if (IRKE_COMMON::MONITOR >= 4) {
                        string central_kmer_string = decode_kmer_from_intval(central_kmer, kmer_length-2);
                        cerr << "input kmer: " << seq << ", central: " << central_kmer_string << endl;
                        
                        string right_central_kmer = decode_kmer_from_intval(get_central_right_kmer(kmer, kmer_length), kmer_length-2);
                        string left_central_kmer = decode_kmer_from_intval(get_central_left_kmer(kmer, kmer_length), kmer_length-2);
                        
                        cerr << "left central kmer: " << left_central_kmer << endl;
                        cerr << "right central kmer: " << right_central_kmer << endl << endl;
                        
                    }
                    
                    // partition kmers according to central kmer value and thread number
                    // so all kmers with common core sequence end up on the same thread.
                    // note, by virtue of this, all extensions for a given kmer should be
                    // accessible via the same node. (idea of Bill Long @ Cray)
                    unsigned int node_for_central_kmer = THIS_MPI_NODE_ID;
                    if (NUM_MPI_NODES > 1) {
                        node_for_central_kmer = get_node_for_central_kmer(central_kmer, kmer_length-2);
                        // otherwise don't waste time looking it up.
                    }
                    if (node_for_central_kmer == THIS_MPI_NODE_ID) {
                        
                        // kmer belongs to this thread.
                        if (WRITE_KMER_FILES) {
                            filewriter << '>' << fe.get_header() << endl << seq << endl;
                        }
                        
                        
                        //kcounter.add_kmer(kmer, count);
                        add_kmer_to_kcounter(kcounter, kmer, count, kmerReaderLog);
                        
                    }
                    else {
                        
                        if (MONITOR_MPI_COMMUNICATION) 
                            kmerReaderLog << "SEND:  Node[" << THIS_MPI_NODE_ID << "] sending 'kmer add' request to node: [" << node_for_central_kmer << "] " << endl;
                        
                        
                        // send it over
                        kmer_int_type_t buffer[5] = {0};
                        buffer[0] = ADD_KMER;
                        buffer[1] = kmer;
                        buffer[2] = count;
                        buffer[3] = THIS_MPI_NODE_ID;
                        buffer[4] = node_for_central_kmer;
                        
                        MPI_Send(&buffer,             /* message buffer */
                                 5,                 /* one data item */
                                 MPI_UNSIGNED_LONG_LONG, /* data item is an integer */
                                 node_for_central_kmer, /* to who we send the message to */
                                 GENERAL_WORK_REQUEST,           /* user chosen message tag */
                                 MPI_COMM_WORLD);   /* default communicator */
                        
                        
                        /*  confirmatory call not essential. Used for debugging.

                        // verify receipt before proceeding
                        MPI_Recv(&buffer,             
                                 5,                 
                                 MPI_UNSIGNED_LONG_LONG, 
                                 node_for_central_kmer,
                                 GENERAL_WORK_RESPONSE, 
                                 MPI_COMM_WORLD, 
                                 &MPI_status);   
                        
                        if (MONITOR_MPI_COMMUNICATION)
                            kmerReaderLog << "RECV: Node[" << THIS_MPI_NODE_ID << "] received " << enum_description[buffer[0]] << " from node[" << node_for_central_kmer << "]" << endl;
                        
                        
                        if (buffer[0] != SENDING_KMER_ADD_CONFIRM){
                            cerr << "ERROR: expected SENDING_KMER_ADD_CONFIRM but received " << enum_description[buffer[0]] << endl;
                            
                        }
                        */
                        
                    }
                    
                    
                }
            }


            cerr << "Node[" << THIS_MPI_NODE_ID << "] is Done populating kmers." << endl;
            if (MONITOR_MPI_COMMUNICATION) 
                kmerReaderLog << "Node[" << THIS_MPI_NODE_ID << "] is Done populating kmers." << endl;
            
            THIS_NODE_DONE = true;
            
            
            if (THIS_MPI_NODE_ID == 0) {
                

	      if (NUM_MPI_NODES == 1) {
		// no reason to run kmer server
		cerr << "** Phase 1: Only 1 MPI node, no reason to do MPI communication. Skipping run_MPI_master_all_completion_check())" << endl;
		
	      }
	      else {
                cerr << "Phase 1: Master node running MPI_completion check.\n";
                
                run_MPI_master_all_completion_check(1);
                
	      }
	    }
            
        }
    }
    
    // add barrier here.
    if (THIS_MPI_NODE_ID == 0) 
        cerr << "REACHED BARRIER." << endl;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    time_t kmer_db_build_end_time;
    time(&kmer_db_build_end_time);
    
    
    if (THIS_MPI_NODE_ID == 0) {
        
        cerr << "BEYOND BARRIER, on to next phase. " << endl;
        
        
        time_t time_to_build_kmer_db =  kmer_db_build_end_time - kmer_db_build_start_time;
            
        cerr << "TIMING_MPI_NODES: " << NUM_MPI_NODES << " TIME_TO_BUILD_KMER_DBS: " << time_to_build_kmer_db << " s." << endl;
	timings_text << "TIMING_REPORT: RANKS: " << NUM_MPI_NODES << " THREADS: " << IRKE_COMMON::NUM_THREADS << " TIME_TO_BUILD_KMER_DBS: " << time_to_build_kmer_db;
    }
    
    
    THIS_NODE_DONE = false;
    ALL_NODES_DONE = false;
    
    if (args.isArgSet("--test_MPI")) {
        // debugging
        test_MPI(kcounter);
        cerr << "Done with testing, now MPI_Finalizing()" << endl;
        MPI_Finalize();        
        exit(0);
    }
    

    // error-kmer pruning
    if (prune_error_kmers) {

        
        if (THIS_MPI_NODE_ID == 0) 
            cerr << "Kmer db size before pruning: " << kcounter.size() << endl;

        kcounter.prune_kmer_extensions(min_ratio_non_error);
        
        if (THIS_MPI_NODE_ID == 0) 
            cerr << "Kmer db size after pruning: " << kcounter.size() << endl;
        

        MPI_Barrier(MPI_COMM_WORLD);
    
        time_t pruning_end_time;
        time(&pruning_end_time);
        
        if (THIS_MPI_NODE_ID == 0) {

            time_t time_for_pruning =  pruning_end_time - kmer_db_build_end_time;
            cerr << "TIMING_MPI_NODES: " << NUM_MPI_NODES << " TIME_TO_PRUNE_ERRORS_FROM_KMER_DBS: " << time_for_pruning << " s." << endl;
	    
	    timings_text << " TIME_TO_PRUNE_ERRORS_FROM_KMER_DBS: " << time_for_pruning;

	}


    }
    

    //------------------------------------------------------------------
    //-- Two parallel processes:
    //    1.   do the contig assembly
    //    2.   act as a server to handle requests from other MPI nodes.
    //------------------------------------------------------------------
    

    time_t contig_building_start_time;
    time(&contig_building_start_time);
       

    /* debugging
    omp_set_num_threads(IRKE_COMMON::NUM_THREADS);
    cerr << "setting num threads to: " << omp_get_num_threads() << endl;
    
    int my_thread_id;
#pragma omp parallel for private (my_thread_id) schedule (static, 1)
    for (int i = 1; i <= 50; i++) {
        my_thread_id = omp_get_thread_num();
        cerr << i << ", thread: " << my_thread_id << endl;
    }
    */
    

   
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section 
        {
                        

            if (is_assembly_node) {


                /* debugging.
                omp_set_num_threads(IRKE_COMMON::NUM_THREADS);
                cerr << "setting num threads to: " << omp_get_num_threads() << endl;
            
                int my_thread_id;
#pragma omp parallel for private (my_thread_id) schedule (static, 1)
                for (int i = 1; i <= 50; i++) {
                    my_thread_id = omp_get_thread_num();
                    cerr << i << ", thread: " << my_thread_id << endl;
                }
                */

                string contig_outfilename = get_MPI_node_filename(THIS_MPI_NODE_ID);
                ofstream contig_writer;
                contig_writer.open(contig_outfilename.c_str());
    
                cerr << "Writing contigs to: " << contig_outfilename << endl;
                
                vector<Kmer_Occurence_Pair> kmers = kcounter.get_kmers_sort_descending_counts(); // not actually sorted!!! //FIXME: hacky use of global sort flag.
                

                // perform assembly using multithreading
                omp_set_num_threads(IRKE_COMMON::NUM_THREADS);
                cerr << "Contig assembly: set OMP threads to: " << IRKE_COMMON::NUM_THREADS << ", and num OMP threads current setting is: " << omp_get_num_threads() << endl;
                int myTid;
                
                map<unsigned int, bool> thread_seen;


                #pragma omp parallel for private (myTid) schedule (dynamic)
                for (unsigned int j = 0; j < kmers.size(); j++) {

                    myTid = omp_get_thread_num();
                    if (thread_seen.find(myTid) == thread_seen.end()) {
                        #pragma omp critical
                        {
                            cerr << "-got contig building thread: " << myTid << endl;
                            thread_seen[ myTid ] = true;
                        }
                    }
                    
                    bool local_debug = false;
                    
                    
                    if (local_debug) {

                        if (THIS_MPI_NODE_ID == 0) 
                            cerr << "thread: " << myTid << " examining kmer: " << j << endl;
                    }
                    
                    kmer_int_type_t kmer = kmers[j].first;
                    unsigned int kmer_count = kcounter.get_kmer_count(kmer);
                    
                    if (! is_good_seed_kmer(kcounter, kmer, kmer_count, kmer_length, MIN_CONNECTIVITY_RATIO)) {
                        continue;
                    }
                    
                    // build draft contig.
                    unsigned int total_counts;
                    vector<kmer_int_type_t> joined_path = build_inchworm_contig_from_seed(kmer, kcounter, MIN_CONNECTIVITY_RATIO, total_counts);
                    

                    // now use this draft contig to select a new seed:
                    kmer_int_type_t new_seed = extract_best_seed(joined_path, kcounter, MIN_CONNECTIVITY_RATIO);
                    
                    if (new_seed == 0) {
                        continue; // must have been zapped by another thread
                    }
            
                    // nicely polished new inchworm contig.
                    joined_path = build_inchworm_contig_from_seed(new_seed, kcounter, MIN_CONNECTIVITY_RATIO, total_counts);
                    
                    string sequence = reconstruct_path_sequence(kcounter, joined_path);
                                        
                    unsigned int contig_length = sequence.length();
                    
                    if (contig_length >= MIN_ASSEMBLY_LENGTH) { // && avg_cov >= MIN_ASSEMBLY_COVERAGE) {
                        
		        #pragma omp critical
                        contig_writer << sequence << endl;     
                        
                    }
                    
                    // remove the kmers.
                    zap_kmers(kcounter, joined_path);
                    
                    
                    /*
                      
                    kmer_int_type_t new_seed = extract_best_seed(joined_path, kcounter, min_connectivity);
                    
                    if (new_seed == 0) {
                    continue; // must have been zapped by another thread
                    }
                    
                    
                    
                    */
                    
                    
                }
                contig_writer.close();
            } // end of assembly-node only

            // all nodes run this.
            THIS_NODE_DONE = true;
            cerr << "NODE: " << THIS_MPI_NODE_ID << " IS DONE ASSEMBLING." << endl << endl;
            
            if (THIS_MPI_NODE_ID == 0) {

	      if (NUM_MPI_NODES == 1) {
		// no reason to run kmer server
		cerr << "** Phase 2: Only 1 MPI node, no reason to do MPI communication. exiting run_MPI_kmer_server()" << endl;
	      }
	      else {
                // master checks other nodes to see if they're done.
                run_MPI_master_all_completion_check(2);
	      }
            }
            
        } // end of assembly OMP section.
        
        #pragma omp section
        {
            //---------------------------
            // Serve the other MPI nodes.
            //---------------------------
            
            //if (is_kmer_server_node) {

	  if (NUM_MPI_NODES == 1) {
	    // no reason to run kmer server
	    cerr << "** Phase 2: Only 1 MPI node, no reason to do MPI communication. Skipping run_MPI_kmer_server()" << endl;
	    
	  }
	  else {
                run_MPI_kmer_server(kcounter);
                // this runs until ALL_NODES_DONE is true.
                // master process sets this 'all done' for all slaves in above under run_MPI_master_all_completion_check()
             //}
	  }
        }
    } // end of omp sections.
    

    time_t contig_building_end_time;
    time(&contig_building_end_time);
    
    time_t time_to_build_contigs = contig_building_end_time - contig_building_start_time;
    if (THIS_MPI_NODE_ID == 0) {
        cerr << "TIMING_MPI_NODES: " << NUM_MPI_NODES << " TIME_TO_BUILD_CONTIGS: " << time_to_build_contigs << " s." << endl;
    
	timings_text << " TIME_TO_BUILD_CONTIGS: " << time_to_build_contigs;
    }
    

    //-------------------------------------------------------
    //-- Pull together final non-redundant inchworm assembly
    //-------------------------------------------------------

    // only master does it.
    if (THIS_MPI_NODE_ID == 0) {


        time_t harvesting_contigs_start_time;
        time(&harvesting_contigs_start_time);
        

        map<unsigned long long, bool> seen_contig_already;
        
        unsigned int INCHWORM_ASSEMBLY_COUNTER = 0;

        int assembly_start_node = 0;
        int assembly_end_node = NUM_MPI_NODES - 1;
        
        if (MASTER_SLAVE_MODE) {
            assembly_end_node = 0;
        }
        else if (NUM_MPI_CONTIG_BUILDERS > 0) {
            assembly_end_node = NUM_MPI_CONTIG_BUILDERS -1;
        }
        
        
        for (unsigned int i=assembly_start_node; i <= assembly_end_node; i++) {
            
            string tmp_contig_file = get_MPI_node_filename(i);
            
            string sequence;
            ifstream tmpreader(tmp_contig_file.c_str());;
            
            while (! tmpreader.eof()) {
                
                tmpreader >> sequence;


                if (tmpreader.eof()) // apparently only happens on the read after the last line is read.
                    break;
                
                
                unsigned int contig_length = sequence.length();
                unsigned int contig_hash = generateHash(sequence);
                
                if (! seen_contig_already[contig_hash]
                    ) {
                    
                    seen_contig_already[contig_hash] = true;
                    
                    INCHWORM_ASSEMBLY_COUNTER++;
                    
                    stringstream headerstream;
                    
                    
                    headerstream << ">a" << INCHWORM_ASSEMBLY_COUNTER 
                        //<< ";" << avg_cov 
                        //       << " total_counts: " <<  total_counts 
                        //<< " Fpath: " << selected_path_n_pair_forward.second << " Rpath: " << selected_path_n_pair_reverse.second 
                        //         << " Seed: " << kmer_count 
                                 << " K: " << kmer_length
                                 << " length: " << sequence.length();
                    
                    string header = headerstream.str();
                    
                    sequence = add_fasta_seq_line_breaks(sequence, 60);
                    
                    cout << header << endl << sequence << endl;      
                    
                    
                }
                
                
                
            }
            
            tmpreader.close();
            if (! KEEP_TMP_FILES) {
                remove(tmp_contig_file.c_str());
            }
        }


        time_t harvesting_contigs_end_time;
        time(&harvesting_contigs_end_time);
        
        time_t contig_harvesting_time = harvesting_contigs_end_time - harvesting_contigs_start_time;

        cerr << "TIMING_MPI_NODES: " << NUM_MPI_NODES << " TIME_TO_HARVEST_CONTIGS: " << contig_harvesting_time << " s." << endl;
	
	timings_text << " TIME_TO_HARVEST_CONTIGS: " << contig_harvesting_time;

	time_t prog_end_time;
	time(&prog_end_time);

	time_t prog_runtime = prog_end_time - prog_start_time;
	
	timings_text << " TIME_TO_RUN_PROG: " << prog_runtime;

	cerr << "DONE." << endl;
	cerr << endl << endl << endl << timings_text.str() << endl << endl << endl << endl;

    }
    

    // all done. :)
    
    MPI_Finalize();

    return(0);
}


string usage (ArgProcessor args) {
    
    bool show_advanced = args.isArgSet("--show_advanced");
    
    stringstream usage_info;
    
    usage_info 
        << endl << endl
        << "Usage:" << endl
        << "  MPIinchworm --kmers <filename.fasta> [opts] " << endl 

        
        << "**Required" << endl
        << "  --kmers  <str>             " << ":fasta file containing kmers and abundance values (from jellyfish)" << endl
        << "    or" << endl
        << "  --reads <str>              " << ":fasta file containing reads" << endl
        << endl
        << "  --token <str>              " << ":a unique alphanumeric to prevent tmp files from clobbering each other from different distributed processes." << endl
        << endl;
    
  usage_info
      << endl
      << "** General opts" << endl
      << "  -K <int>                         " << ":kmer length (default: 25, meaning 24mer overlaps)  max = 32 (stored as 64-bit integers, 2-base char encoding)" << endl
      //<< "  --minKmerCount <int>             " << ":min kmer count, lower abundant kmers are pruned from graph (default: 1)" << endl
      << "  -L <int>                         " << ":min contig length to be reported (default: 25)" << endl
      //<< "  --min_assembly_coverage  <int>   " << ":min kmer coverage of an assembled contig to be reported (default: 2)" << endl
      << "  --DS                             " << ":double-stranded RNA-Seq mode (not strand-specific)" << endl
      << "  --num_threads <int>              " << ":number of threads to use. (by default, uses OMP_NUM_THREADS env var setting)" << endl
      << "  --num_mpi_contig_builders <int>  " << ":indicate number of mpi contig builder nodes. This partitions nodes into disjoint sets of contig builders and kmer servers, where remaining mpi nodes are kmer servers. (default: -1, all are both contig builders and kmer servers)" << endl
      << "  --no_prune_error_kmers           " << ":disable pruning of kmers that occur at below --min_ratio_non_error " << endl
      << "  --min_ratio_non_error <float>    " << ":min ratio for low/high alternative extension that is not an error (default: 0.05)" << endl
      


      ;
    
  if (show_advanced) {
        usage_info 
            << "    * advanced opts:" << endl
            << "    --min_seed_entropy <float> " << ":min seed kmer entropy (default: 1.5)" << endl
            << "    --min_seed_coverage <int> " << ":min seed kmer coverage (default: 2)" << endl;
  
        usage_info << "   *debugging opts " << endl
                   << "   --test_MPI               :just test MPI calls." << endl
                   << "        --max_test_kmers <int>   : number of kmers to iterate through in MPI call testing. (default: " << MAX_TEST_KMERS << ")" << endl
                   << "  --master_slave_mode       :only the master builds contigs, all other nodes are kmer servers." << endl
                   << "  --write_kmer_files        :each node writes out its subset of stored kmers to a separate file." << endl
                   << "  --monitor_MPI             :verbosely describe MPI communications underway." << endl
                   << "  --keep_tmp_files          :retain temporary files. " << endl
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

kmer_int_type_t get_central_kmer(kmer_int_type_t kmer, unsigned int kmer_length) {

    // given ABCDE, want BCD

    kmer >>=2; // remove last nucleotide
    
    kmer_int_type_t kmer_mask = pow(2,2*( (kmer_length-1) -1) ) -1; // remove first nucleotide of the resulting (kmer-1) to get the core seq
        
    kmer_int_type_t central_kmer = kmer & kmer_mask;
    
    return(central_kmer);
    
}


kmer_int_type_t get_central_right_kmer(kmer_int_type_t kmer, unsigned int kmer_length) {
    
    // given ABCDE, want CDE

    kmer_int_type_t kmer_mask = pow(2,2*(kmer_length-2))-1; // remove first two nucleotides of kmer
        
    kmer_int_type_t central_kmer = kmer & kmer_mask;
    
    return(central_kmer);
    
}

kmer_int_type_t get_central_left_kmer(kmer_int_type_t kmer, unsigned int kmer_length) {

    // given ABCDE, want ABC

    kmer >>= 4; // shift out the last two nucleotides.

    return(kmer);
}


bool is_good_seed_kmer(KmerCounter& kcounter, kmer_int_type_t kmer, unsigned int kmer_count, unsigned int kmer_length, 
                       float min_connectivity) {

    if (kmer_count == 0) {
        return(false);
    }
    
    if (kmer == revcomp_val(kmer, kmer_length)) {
        // palindromic kmer, avoid palindromes as seeds
                
        if (IRKE_COMMON::MONITOR >= 10) {
            cerr << "SEED kmer: " << kcounter.get_kmer_string(kmer) << " is palidnromic.  Skipping. " << endl;
        }
        
        return(false);
    }
    
        
    if (kmer_count < MIN_SEED_COVERAGE) {
        if (IRKE_COMMON::MONITOR >= 10) {
            cerr << "-seed has insufficient coverage, skipping" << endl;
        }
        
        return(false);
    }
    
    float entropy = compute_entropy(kmer, kmer_length);
        
        
    if (entropy < MIN_SEED_ENTROPY) {
        
        if (IRKE_COMMON::MONITOR >= 10) {
            cerr << "-skipping seed due to low entropy: " << entropy << endl;
        }
        
        return(false);
    }
        

    // got this far, so kmer is fine as a seed
    return(true);
}

vector<kmer_int_type_t>  build_inchworm_contig_from_seed(kmer_int_type_t kmer, KmerCounter& kcounter, float min_connectivity, unsigned int& total_counts) {

    unsigned int kmer_count = kcounter.get_kmer_count(kmer); 
        
    unsigned int kmer_length = kcounter.get_kmer_length();
    
    // track those kmers included in growing path.
    Kmer_visitor visitor(kmer_length, DOUBLE_STRANDED_MODE);
    
    //-----------------------
    /* Extend to the right */
    //-----------------------
    
    vector<kmer_int_type_t> forward_path = inchworm(kcounter, 'F', kmer, visitor, min_connectivity); 
    
   
    if (IRKE_COMMON::MONITOR >= 2) {
        cerr << "Forward path contains: " << forward_path.size() << " kmers. " << endl;
    
        for (unsigned int i = 0; i < forward_path.size(); i++) {
            kmer_int_type_t kmer = forward_path[i];
            cerr << "\tForward path kmer: " << kcounter.get_kmer_string(kmer) << endl;
        }
    }
        
    //----------------------
    /* Extend to the left */ 
    //----------------------
    
    vector<kmer_int_type_t> reverse_path = inchworm(kcounter, 'R', kmer, visitor, min_connectivity);
    
    
    if (IRKE_COMMON::MONITOR >= 2) {
        cerr << "Reverse path contains: " << reverse_path.size() << " kmers. " << endl;
        for (unsigned int i = 0; i < reverse_path.size(); i++) {
            cerr  << "\tReverse path kmer: " << kcounter.get_kmer_string(reverse_path[i]) << endl; 
        }
    }
                    
                
    vector<kmer_int_type_t> joined_path = join_forward_n_reverse_paths(reverse_path, kmer, forward_path);

    return(joined_path);

}


vector<kmer_int_type_t> join_forward_n_reverse_paths(vector<kmer_int_type_t>& reverse_path, kmer_int_type_t seed_kmer_val, vector<kmer_int_type_t>& forward_path) {
    
    vector<kmer_int_type_t> joined_path;
        
    // want reverse path in reverse order
        
    for (int i = reverse_path.size()-1; i >= 0; i--) {
        joined_path.push_back( reverse_path[i] );
    }
        
    // add seed kmer
    joined_path.push_back(seed_kmer_val);
        
    // tack on the entire forward path.
        
    for (unsigned int i = 0; i < forward_path.size(); i++) {
        joined_path.push_back( forward_path[i] );
    }
        
    return(joined_path);
}

vector<kmer_int_type_t> inchworm (KmerCounter& kcounter, char direction, kmer_int_type_t kmer, Kmer_visitor& visitor, float min_connectivity) {

    vector<kmer_int_type_t> growing_path;
    
    int kmer_length = kcounter.get_kmer_length();

    while (true) {
        kmer_int_type_t next_kmer_core;
        
        if (direction == 'F') {
            next_kmer_core = get_central_right_kmer(kmer, kmer_length);
        }
        else {
            next_kmer_core = get_central_left_kmer(kmer, kmer_length);
        }
        
        // figure out which node
        unsigned int node_for_central_kmer = get_node_for_central_kmer(next_kmer_core, kmer_length-2);
        
        kmer_int_type_t best_extension = 0;
        
        if (node_for_central_kmer == THIS_MPI_NODE_ID) {
            // can pull the next kmer counts from our own hash table
            vector<Kmer_Occurence_Pair> kmer_candidates;
            if (direction == 'F') {
                // forward search
                kmer_candidates = kcounter.get_forward_kmer_candidates(kmer);
            }
            else {
                // reverse search
                kmer_candidates = kcounter.get_reverse_kmer_candidates(kmer);
            }
            if (kmer_candidates.size()) {
                best_extension = kmer_candidates[0].first;
            }
            
        }
        else {
            // need to send MPI request to the node that has that kmer count set.
            
            best_extension = get_greedy_extension_from_MPI_node(node_for_central_kmer, kmer, direction, kmer_length);
            
        }
        
        if (best_extension == 0) {
            break;
        }
        else if (visitor.exists(best_extension)) {
            break;
        }
        else {
            visitor.add(best_extension);
            growing_path.push_back(best_extension);
            kmer = best_extension;
        }
    }
    
    return(growing_path);
}


unsigned int get_node_for_central_kmer (kmer_int_type_t central_kmer, unsigned int kmer_length) {
    

    kmer_int_type_t canonical_central_kmer = central_kmer;
    kmer_int_type_t rev_central_kmer = revcomp_val(central_kmer, kmer_length);
    
    if (rev_central_kmer < canonical_central_kmer) {
        canonical_central_kmer = rev_central_kmer;
    }
    

    unsigned int node_for_kmer = canonical_central_kmer % NUM_MPI_NODES;
    
    if (IRKE_COMMON::MONITOR >= 4) {
        cerr << "Kmer: " << decode_kmer_from_intval(central_kmer, kmer_length) << " or " 
             << decode_kmer_from_intval(canonical_central_kmer, kmer_length) << " assigned to node: " << node_for_kmer << endl; 
    }
    
    // all nodes are kmer servers
    return(node_for_kmer);
    
}


string reconstruct_path_sequence(KmerCounter& kcounter, vector<kmer_int_type_t>& path) {
        
    if (path.size() == 0) {
        return("");
    }
        
    string seq = kcounter.get_kmer_string(path[0]);
    //cov_counter.push_back( kcounter.get_kmer_count(path[0]) );
        
    for (unsigned int i = 1; i < path.size(); i++) {
        string kmer = kcounter.get_kmer_string(path[i]);
        seq += kmer.substr(kmer.length()-1, 1);
                
        //cov_counter.push_back( kcounter.get_kmer_count(path[i]) );
    }
    
    return(seq);
}



void zap_kmers(KmerCounter& kcounter, vector<kmer_int_type_t>& kmer_path) {
    
    unsigned int kmer_length = kcounter.get_kmer_length();

    kmer_int_type_t buffer[5] = {0};
    buffer[0] = CLEAR_KMER;

    // exclude kmers built into current contig.
    for (unsigned int i = 0; i < kmer_path.size(); i++) {
        
        kmer_int_type_t kmer = kmer_path[i];

        // determine MPI node placement
        
        kmer_int_type_t central_kmer = get_central_kmer(kmer, kmer_length);
        
        unsigned int MPI_node_id = get_node_for_central_kmer(central_kmer, kmer_length-2);
        
        if (MPI_node_id == THIS_MPI_NODE_ID) {
        
            kcounter.clear_kmer(kmer);
        }
        else {
            // send message to that node to zap the kmer
           
            buffer[1] = kmer;
            buffer[3] = THIS_MPI_NODE_ID; // from
            buffer[4] = MPI_node_id; // to
            
            MPI_Send(&buffer,             /* message buffer */
                     5,                 /* one data item */
                     MPI_UNSIGNED_LONG_LONG, /* data item is an integer */
                     MPI_node_id, /* to who we just received from */
                     GENERAL_WORK_REQUEST,           /* user chosen message tag */
                     MPI_COMM_WORLD);   /* default communicator */
            
        }
        
    }
}

    
void run_MPI_master_all_completion_check(int phase_val) {

    // if not the master, just return...  
    if (THIS_MPI_NODE_ID != 0) {
        return;
    }
    

    cerr << "-running MPI master all completion check" << endl;

    ofstream filewriter;
    if (MONITOR_MPI_COMMUNICATION) {
        stringstream filenamebuilder;
        filenamebuilder << "master.completion_check." << phase_val << ".log";
        filewriter.open(filenamebuilder.str().c_str());
    }


    map<unsigned int,bool> unfinished_nodes;
    map<unsigned int,bool>::iterator it;
    
    for (int i = 1; i < NUM_MPI_NODES; i++) { // do not include self (master) node - already know master is done.
        unfinished_nodes[i] = true;
    }
    

    kmer_int_type_t buffer[5] = {0};
               
    MPI_Status MPI_status; 
    
    while (! unfinished_nodes.empty()) {
        
        vector<unsigned int> finished_nodes;
        
        for (it = unfinished_nodes.begin(); it != unfinished_nodes.end(); it++) {
            
            unsigned int node_id = it->first;
            // check to see if it's finished.
                       
            buffer[0] = GET_FINISHED_STATUS;
            buffer[1] = 0; // init to false
            
            buffer[3] = node_id; // from
            buffer[4] = THIS_MPI_NODE_ID; // to

            
            cerr << "asking node: " << node_id << " if finished... ";
            if (MONITOR_MPI_COMMUNICATION) {
                filewriter << ">> asking GET_FINISHED_STATUS from node[" << node_id << "] " << endl;
            }


            // require MPI call here to that node to ask if its done yet.
            MPI_Send(&buffer,             /* message buffer */
                     5,                 /* one data item */
                     MPI_UNSIGNED_LONG_LONG, /* data item is an integer */
                     node_id, /* to who we just received from */
                     GENERAL_WORK_REQUEST,           /* user chosen message tag */
                     MPI_COMM_WORLD);   /* default communicator */
            
            MPI_Recv(&buffer,             /* message buffer */
                     5,                 /* one data item */
                     MPI_UNSIGNED_LONG_LONG, /* data item is an integer */
                     node_id, /* to who we just received from */
                     GENERAL_WORK_RESPONSE,           /* user chosen message tag */
                     MPI_COMM_WORLD, /* default communicator */
                     &MPI_status);   

            cerr << "?? IS_FINISHED? from master to node[" << node_id << "] is: " << buffer[1] << ", response msg: " << enum_description[buffer[0]] 
                 << ", audit_from: " << buffer[3]
                 << ", audit to: " << buffer[4] << endl;

            if (MONITOR_MPI_COMMUNICATION) {
                filewriter << "<< received GET_FINISHED_STATUS from node[" << node_id << "] = " << enum_description[buffer[0]] << endl;
            }
            
            if (buffer[0] != SENDING_FINISHED_STATUS || buffer[0] == ERROR_ENCOUNTERED || buffer[3] != node_id || buffer[4] != THIS_MPI_NODE_ID) {
                cerr << endl << "** ERROR ** getting finished status from node: " << node_id << endl << endl;
            }
            else {
                if (buffer[1] == 1) {
                    // node is finished.
                    finished_nodes.push_back(node_id);
                }
            }
            cerr << "\tchecking node: " << node_id << " if finished: " << buffer[1] << endl;
            
        }
        
        for (int i =0; i < finished_nodes.size(); i++) {
            unfinished_nodes.erase(finished_nodes[i]);
        }
        if (! unfinished_nodes.empty()) {
            cerr << "resting between polls for " << MPI_SLEEP_INTERVAL << endl;
            sleep(MPI_SLEEP_INTERVAL);
            
        }
        
        cerr << "Number of unfinished nodes is currently: " << unfinished_nodes.size() << endl;
    }
    
    // send any waiting nodes the 'all done' signal.
    buffer[0] = SET_ALL_NODES_DONE;
    
    int first_node = 0;
    if (MASTER_SLAVE_MODE)
        first_node = 0; // master is not running a kmer server, so no need to send a msg to it.  probably doesnt matter though. // yes it is!!!
         
    for (int i = first_node; i < NUM_MPI_NODES; i++) { // include msg to master itself since the other thread is likely waiting for a msg.
    
        cerr << "Signalling all done to node: " << i << endl;
        buffer[3] = THIS_MPI_NODE_ID; // from
        buffer[4] = i; // to
        
        MPI_Send(&buffer,             /* message buffer */
                 5,                 /* one data item */
                 MPI_UNSIGNED_LONG_LONG, /* data item is an integer */
                 i, /* to who we just received from */
                 GENERAL_WORK_REQUEST,           /* user chosen message tag */
                 MPI_COMM_WORLD);   /* default communicator */
    }
    
    
}

void test_MPI(KmerCounter& kcounter) {

    // master (rank=0) loops through local kmer set and queries the other nodes for extensions.

    unsigned int kmer_length = kcounter.get_kmer_length();
    
    bool is_assembly_node = true;
    bool is_kmer_server_node = true;
    
    if (MASTER_SLAVE_MODE) {
        if (THIS_MPI_NODE_ID == 0) {
            is_assembly_node = true;
            is_kmer_server_node = false;
        }
        else {
            is_assembly_node = false;
            is_kmer_server_node = true;
        }
    }
    

    
    #pragma omp parallel sections num_threads(2)
    {
         #pragma omp section 
        {
            if (is_assembly_node) {
                vector<Kmer_Occurence_Pair> kmers = kcounter.get_kmers_sort_descending_counts();
                
                int counter = 0;
                for (unsigned int i = 0; i < kmers.size(); i++) {
                    
                    kmer_int_type_t kmer = kmers[i].first;
                    
                    if (THIS_MPI_NODE_ID == 0 && IRKE_COMMON::MONITOR >= 4) {
                        cerr << "Kmer is: " << decode_kmer_from_intval(kmer, kmer_length) << " with count: " << kcounter.get_kmer_count(kmer) << endl;
                    }
                    

                    kmer_int_type_t next_kmer_core = get_central_right_kmer(kmer, kmer_length);
                    unsigned int node_for_central_kmer = get_node_for_central_kmer(next_kmer_core, kmer_length-2); 
                    
                    if (node_for_central_kmer == THIS_MPI_NODE_ID) {
                        // only asking other nodes for kmers.
                        continue;
                    }
                    
                    // need to send MPI request to the node that has that kmer count set.
                    
                    kmer_int_type_t greedy_right_kmer = get_greedy_extension_from_MPI_node(node_for_central_kmer, kmer, 'F', kmer_length);

                    if (THIS_MPI_NODE_ID == 0) 
                        cerr << "\t\tRIGHT kmer stored at node: " << node_for_central_kmer << " is " 
                             << decode_kmer_from_intval(greedy_right_kmer, kmer_length) << endl;
                    
                    
                    kmer_int_type_t greedy_left_kmer = get_greedy_extension_from_MPI_node(node_for_central_kmer, kmer, 'R', kmer_length);
                    
                    if (THIS_MPI_NODE_ID == 0) 
                        cerr << "\t\tLEFT kmer stored at node: " << node_for_central_kmer << " is " 
                             << decode_kmer_from_intval(greedy_left_kmer, kmer_length) << endl;
                    
                    
                    counter++;
                    
                    if (counter > MAX_TEST_KMERS) { 
                        
                        if (THIS_MPI_NODE_ID == 0)
                            cerr << "** Stopping now, reached MAX_TEST_KMERS, counter at: " << counter << endl;
                        
                        break; 
                    } 
                    
                }
                
                cerr << "NODE: [" << THIS_MPI_NODE_ID << "] is done with kmer extension test." << endl;
                
                if (THIS_MPI_NODE_ID == 0) { // master
                    run_MPI_master_all_completion_check(3);
                }
                
                
            }
            
            THIS_NODE_DONE = true;
        }
        
        #pragma omp section 
        {
            
            if (is_kmer_server_node) {
                
                run_MPI_kmer_server(kcounter);
                
            }
            
        }   
    }
}



void run_MPI_kmer_server(KmerCounter& kcounter) {
  
  KMER_SERVER_PHASE++;

  cerr << "** Launching run_MPI_kmer_server() for node: " << THIS_MPI_NODE_ID << ", kmer_server_phase: " << KMER_SERVER_PHASE << endl;
  
    stringstream filenamebuilder;
    filenamebuilder << "mpi_kmer_server.mpi_" << THIS_MPI_NODE_ID << ".phase" << KMER_SERVER_PHASE << ".log";
    
    ofstream kmerServerLog;
    
    if (MONITOR_MPI_COMMUNICATION) {
        kmerServerLog.open(filenamebuilder.str().c_str());
    }   
    
    MPI_Status MPI_status;
    
    while (! ALL_NODES_DONE) {
        
        unsigned long long request[5] = {0};
        
        if (MONITOR_MPI_COMMUNICATION)  {
            kmerServerLog  << "|| MPISERVER_node:[" << THIS_MPI_NODE_ID << "] Status update: serving other mpi nodes." << endl;
        }
        
        /* Receive request */
        
        MPI_Recv(&request,           /* message buffer */
                 5,                 /* one data item */
                 MPI_UNSIGNED_LONG_LONG,        /* of type double real */
                 MPI_ANY_SOURCE,    /* receive from any sender */
                 GENERAL_WORK_REQUEST,       /* any type of message */
                 MPI_COMM_WORLD,    /* default communicator */
                 &MPI_status);          /* info about the received message */
        
        if (MONITOR_MPI_COMMUNICATION)  {
            kmerServerLog << "<< MPISERVER_node:[" << THIS_MPI_NODE_ID << "] received GENERAL_WORK_REQUEST (" << enum_description[request[0]] << "), value (" << request[1] << "), from node: [" << MPI_status.MPI_SOURCE << "]" << ", audit_from: " << request[3] << ", audit to: " << request[4] << endl;
        }
        
        bool response_required = false;
        
        unsigned long long response[5] = {0};
        
        
        
        if (request[0] == GET_RIGHT_GREEDY_EXTENSION) {
            response_required = true;
            response[0] = SENDING_RIGHT_GREEDY_EXTENSION;
            kmer_int_type_t input_kmer = request[1];
            vector<Kmer_Occurence_Pair> right_extensions_sorted = kcounter.get_forward_kmer_candidates(input_kmer);
            if (right_extensions_sorted.size()) {
                Kmer_Occurence_Pair extend_kmer_pair = right_extensions_sorted[0];
                kmer_int_type_t extend_kmer = extend_kmer_pair.first;
                unsigned int count = extend_kmer_pair.second;
                response[1] = extend_kmer;
                response[2] = count;
            }
            
        }
        else if (request[0] == GET_LEFT_GREEDY_EXTENSION) {
            response_required = true;
            response[0] = SENDING_LEFT_GREEDY_EXTENSION;
            kmer_int_type_t input_kmer = request[1];
            vector<Kmer_Occurence_Pair> left_extensions_sorted = kcounter.get_reverse_kmer_candidates(input_kmer);
            if (left_extensions_sorted.size()) {
                Kmer_Occurence_Pair extend_kmer_pair = left_extensions_sorted[0];
                kmer_int_type_t extend_kmer = extend_kmer_pair.first;
                unsigned int count = extend_kmer_pair.second;
                response[1] = extend_kmer;
                response[2] = count;
            }
        }
        else if (request[0] == SET_ALL_NODES_DONE) {
            ALL_NODES_DONE = true;
            
        }
        else if (request[0] == GET_FINISHED_STATUS) {
            response_required = true;
            response[0] = SENDING_FINISHED_STATUS;
            response[1] = THIS_NODE_DONE;
        }
        else if (request[0] == CLEAR_KMER) {
            kmer_int_type_t kmer = request[1];
            kcounter.clear_kmer(kmer);
            
        }
        else if (request[0] == ADD_KMER) {
            
            kmer_int_type_t kmer = request[1];
            unsigned int count = (unsigned int) request[2];
            //kcounter.add_kmer(kmer, count);
            add_kmer_to_kcounter(kcounter, kmer, count, kmerServerLog);
            
            //response_required = true;      
            //response[0] = SENDING_KMER_ADD_CONFIRM;
            
            response_required = false;      
            
            
        }
        
        else {
            // not sure what the request is
            response_required = true;
            response[0] = ERROR_ENCOUNTERED;
        }
        
        if (response_required && ! ALL_NODES_DONE) {
            
            // send response
            
            if (MONITOR_MPI_COMMUNICATION)  {
                kmerServerLog << "....>> MPISERVER (sending) node:[" << THIS_MPI_NODE_ID << "] sending response " << enum_description[response[0]] << " to node: [" << MPI_status.MPI_SOURCE << "]" << endl;
                
            }
            
            response[3] = THIS_MPI_NODE_ID; // from 
            response[4] = MPI_status.MPI_SOURCE; // to
            
            MPI_Send(&response,             /* message buffer */
                     5,                 /* one data item */
                     MPI_UNSIGNED_LONG_LONG, /* data item is an integer */
                     MPI_status.MPI_SOURCE, /* to who we just received from */
                     GENERAL_WORK_RESPONSE,           /* user chosen message tag */
                     MPI_COMM_WORLD);   /* default communicator */
            
            if (MONITOR_MPI_COMMUNICATION) {
                kmerServerLog << ">> MPISERVER (SENT) node:[" << THIS_MPI_NODE_ID << "] msg GENERAL_WORK_RESPONSE of " << enum_description[response[0]] 
                              << " value: (" << response[1] << ") sent to node[" << MPI_status.MPI_SOURCE << "]" << endl; 
            }
        }
        
        
    }
    
    cerr << "** Shutting down kmer server at node[" << THIS_MPI_NODE_ID << "] from:  run_MPI_kmer_server()" << endl;
    if (MONITOR_MPI_COMMUNICATION) {
      kmerServerLog << "** Shutting down kmer server at node[" << THIS_MPI_NODE_ID << "] from:  run_MPI_kmer_server()" << endl;
      kmerServerLog.close();
    }
}



kmer_int_type_t get_greedy_extension_from_MPI_node(unsigned int node_for_central_kmer, kmer_int_type_t kmer, char direction, unsigned int kmer_length) {


    if (IRKE_COMMON::MONITOR >= 4) 
        cerr << "** SEND:  Node " << THIS_MPI_NODE_ID << " requests extension from node: " << node_for_central_kmer << " kmer: " << kmer << " direction: " << direction << endl;

    unsigned long long buffer[5] = {0};

    unsigned int WORK_TAG;
    unsigned int EXPECTED_RESPONSE_TAG;

    if (direction == 'F') {
        WORK_TAG = GET_RIGHT_GREEDY_EXTENSION;
        EXPECTED_RESPONSE_TAG = SENDING_RIGHT_GREEDY_EXTENSION;
    }
    else if (direction == 'R') {
        WORK_TAG = GET_LEFT_GREEDY_EXTENSION;
        EXPECTED_RESPONSE_TAG = SENDING_LEFT_GREEDY_EXTENSION;
    }
    else {
        throw(stacktrace() + "\n\nOnly recognize direction of 'F' or 'R'");
    }

    buffer[0] = WORK_TAG;
    buffer[1] = kmer;

    buffer[3] = THIS_MPI_NODE_ID; // from
    buffer[4] = node_for_central_kmer; // to
        
    #pragma omp critical (kmer_extend_request)
    {
        MPI_Send(&buffer,             /* message buffer */
                 5,                 /* one data item */
                 MPI_UNSIGNED_LONG_LONG, /* data item is an integer */
                 node_for_central_kmer, /* to who we just received from */
                 GENERAL_WORK_REQUEST,           /* user chosen message tag */
                 MPI_COMM_WORLD);   /* default communicator */
        
        
        
        /* Receive request */
        
        if (MONITOR_MPI_COMMUNICATION) 
            cerr << "** RECV(waiting): Node [" << THIS_MPI_NODE_ID << "] waiting on extension from node: [" << node_for_central_kmer << "] and looking for response tag: " << enum_description[EXPECTED_RESPONSE_TAG] << endl;
        
        
        MPI_Status MPI_status;
        
        MPI_Recv(&buffer,
                 5,      
                 MPI_UNSIGNED_LONG_LONG, 
                 node_for_central_kmer,
                 GENERAL_WORK_RESPONSE, 
                 MPI_COMM_WORLD,  
                 &MPI_status);    
    }
    
    if (MONITOR_MPI_COMMUNICATION) 
        cerr << "** RECV: Node [" << THIS_MPI_NODE_ID << "] received response tag (" << enum_description[buffer[0]] << ") from node: [" << node_for_central_kmer << "] kmer: " 
             << decode_kmer_from_intval(kmer, kmer_length) << " direction: " << direction << " as " << decode_kmer_from_intval(buffer[1], kmer_length)
             << ", audit_from: " << buffer[3] << ", audit_to: " << buffer[4] << endl;
    


    if (buffer[0] != EXPECTED_RESPONSE_TAG || buffer[0] == ERROR_ENCOUNTERED) {
        // handle error
        cerr << stacktrace() << "ERROR ENCOUNTERED IN RETRIEVING EXTENSION VIA MPI: " 
             << "node[" << THIS_MPI_NODE_ID << "] asked node [" << node_for_central_kmer << "]" << endl
             << "\tfor " << direction << " extension of " << decode_kmer_from_intval(kmer, kmer_length) << endl
             << "\tand received status: " << buffer[0] << " and answer " << decode_kmer_from_intval(buffer[1], kmer_length) << endl 
             << "\tand expected response tag was: " << EXPECTED_RESPONSE_TAG << endl
             << endl;
        
        return(0);
    }
    else {
        return(buffer[1]);
    }
}


string get_MPI_node_filename(unsigned int node_id) {
     
    stringstream filename_builder;
    filename_builder << "tmp." << TOKEN << ".iworm_mpi_node_" << node_id << ".contigs.txt";
    
    return(filename_builder.str());
}

string add_fasta_seq_line_breaks(string& sequence, int interval) {
    
    stringstream fasta_seq;
    
    int counter = 0;
    for (string::iterator it = sequence.begin(); it != sequence.end(); it++) {
        counter++;
        
        fasta_seq << *it;
        if (counter % interval == 0 && (it + 1) != sequence.end()) {
            fasta_seq << endl;
        }
    }
    
    return(fasta_seq.str());
}


kmer_int_type_t extract_best_seed(vector<kmer_int_type_t>& kmer_vec, KmerCounter& kcounter, float min_connectivity) {

    unsigned int kmer_length = kcounter.get_kmer_length();
    
    unsigned int best_kmer_count = 0;
    kmer_int_type_t best_seed;
    
    for (unsigned int i = 0; i < kmer_vec.size(); i++) {
        
        kmer_int_type_t kmer = kmer_vec[i];
        unsigned int count = kcounter.get_kmer_count(kmer);
        
        if (count > best_kmer_count && is_good_seed_kmer(kcounter, kmer, count, kmer_length, min_connectivity)) {
            best_kmer_count = count;
            best_seed = kmer;
        }
    }
    
    if (IRKE_COMMON::MONITOR >= 2) {
        cerr << "Parallel method found better seed: " << kcounter.get_kmer_string(best_seed) << " with count: " << best_kmer_count << endl;
    }

    return(best_seed);
}


void add_kmer_to_kcounter(KmerCounter& kcounter, kmer_int_type_t& kmer, unsigned int& count, ofstream& logstream) {

  /*
  if (MONITOR_MPI_COMMUNICATION) {
    logstream << "NODE[" << THIS_MPI_NODE_ID << "] adding kmer: " << kmer << " with count " << count << " to KmerCounter" << endl;
  }
  */  

  #pragma omp critical (KmerAdding)
  {
    kcounter.add_kmer(kmer, count);

  }

}
