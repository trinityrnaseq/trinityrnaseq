
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "analysis/TranscriptomeGraph.h"
#include "analysis/DNAVector.h"
#include "analysis/CompMgr.h"

/*
  
  #############################
  # overall Chrysalis workflow:
  #############################

  -run bowtie for scaffolding

  -run graphFromFasta
  
  -create the bundled.fasta file
  
  -run reads-to-transcripts
  
  -sort the reads-to-transcripts file.
  
  -write the deBruijn graphs

 */



int KMER_SIZE = 24;
typedef vector<string> string_vec;
static bool NO_CLEANUP = false;

static bool __DEBUG_WELD_ALL = false;

string sort_exec = "sort";

void Execute(const char * command) {
    int ret = system(command);
    if (ret != 0) {
        cout << "COMMAND: " << command << endl;
        cout << "Died with exit code " << ret << endl;
        cout << "Exiting." << endl;
        exit(-1);
    }
    
}

bool Exists(const string & s) 
{
    FILE * p = fopen(s.c_str(), "r");
    if (p != NULL) {
        fclose(p);
        return true;
    }
    // cout << "FATAL ERROR: Could not open file for read: " << s << endl;
    // cout << "Please make sure to enter the correct file name(s). Exiting now." << endl;
    
    return false;
}



void write_deBruijn_graphs(string iwormBundleFasta, string execDir, bool sStrand, int nCPU) {

   string cpp_graph_writer = execDir + "/../Inchworm/bin/FastaToDeBruijn";

   
   string graph_output_filename = iwormBundleFasta + ".deBruijn";
   
   stringstream cpp_graph_cmd;
   cpp_graph_cmd << cpp_graph_writer << "  --fasta " << iwormBundleFasta
       // << "," << reads_file // use both reads and iworm contigs for building graph.
                 << " -K " << KMER_SIZE << " --graph_per_record --threads " << nCPU;
   if (sStrand) {
       cpp_graph_cmd << " --SS ";
   }
   cpp_graph_cmd << " > " << graph_output_filename;
   
   // cerr << "running CMD: [" << i << "] " << cpp_graph_cmd.str() << endl;
   
   Execute(cpp_graph_cmd.str().c_str());
   
}



void write_deBruijn_graphs(vector<string_vec>& bundled, vector<int>& component_ids, ComponentFileMgr& mgr, string outDir, string execDir, bool sStrand, map<int,bool>& pursue_component) {

    string cpp_graph_writer = execDir + "/../Inchworm/bin/FastaToDeBruijn";
    
    int num_finished = 0;
    int num_components = bundled.size();

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < bundled.size(); i++) {

        int component_id = component_ids[i];
        
        if (pursue_component.find(component_id) == pursue_component.end()) {
            // not pursuing it.
            continue;
        }


        string_vec bundled_iworms = bundled[i];
        
        string component_file_basename = mgr.GetFileName(component_id, "");
        
        
        string graph_filename = component_file_basename + ".raw.graph";
        string comp_bundle_file = component_file_basename + ".iworm_bundle";
        
        if (Exists(graph_filename) && ! Exists(comp_bundle_file)) {
            cerr << "-already processed graph for: " << graph_filename << ", must be resuming operations." << endl;
            continue;
        }
        

        // use reads instead if iworm bundles for building de bruijn graph:
        string reads_file = component_file_basename + ".raw.fasta";
        
        stringstream cpp_graph_cmd;
        cpp_graph_cmd << cpp_graph_writer << "  --fasta " << comp_bundle_file 
                   // << "," << reads_file // use both reads and iworm contigs for building graph.
                      << " -K " << KMER_SIZE << " -C " << component_id;
        if (sStrand) {
           cpp_graph_cmd << " --SS ";
        }
        cpp_graph_cmd << " > " << graph_filename;
        
        // cerr << "running CMD: [" << i << "] " << cpp_graph_cmd.str() << endl;
        
        Execute(cpp_graph_cmd.str().c_str());
        
        // cleanup
        if (! NO_CLEANUP) {
            // bundle file served its purpose.
            string cleanup_cmd = "rm " + comp_bundle_file;
            Execute(cleanup_cmd.c_str());
        }
        
        
        #pragma omp atomic
        num_finished++;
        
        #pragma omp critical 
        {

            cerr << "\r" << num_finished << "/" << num_components << " = " << (float)num_finished/num_components*100 << "% of de Bruijn graphs constructed   ";
            //cerr << "Executing: " << graph_cmd.str() << endl;
        }

    }
    
    
    cerr << "-done writing de Bruijn graphs." << endl;
    
    return;
}



void write_iworm_bundle (string filename, vector<string_vec>& bundled, vector<int>& component_ids, ComponentFileMgr& mgr) {

    ofstream ofh;
    ofh.open(filename.c_str());
    
    //ofstream bundle_listing_ofh;
    //bundle_listing_ofh.open("iworm_bundle_file_listing.txt");


    for (int i = 0; i < bundled.size(); i++) {

        int component_id = component_ids[i];
        string_vec bundled_iworms = bundled[i];

        stringstream s;

        s << ">s_" << component_id << endl;
        
        int num_iworms = bundled_iworms.size();
        
        for (int j = 0; j < num_iworms; j++) {
            s << bundled_iworms[j];
            if (j < num_iworms-1) {
                s << "X";
            }
        }
        s << endl;
        
        ofh << s.str();
        
        /*  minimize number of files to be generated.

        // write a local version of it according to the component identifier
        string comp_bundle_file = mgr.GetFileName(component_id, ".iworm_bundle");
        ofstream bundle_ofh;
        bundle_ofh.open(comp_bundle_file.c_str());
        bundle_ofh << s.str();
        bundle_ofh.close();

        bundle_listing_ofh << comp_bundle_file.c_str() << endl;

        */
        
    }

    ofh.close();

    //bundle_listing_ofh.close();
    
    return;
    
    
}

string get_seq_string(DNAVector& d) {

    stringstream s;
    
    for (int i = 0; i < d.isize(); i++) {
        s << d[i];
    }

    return(s.str());
}





int main(int argc,char** argv)
{
    char execPath[512];
    strcpy(execPath, argv[0]);
    execPath[strlen(execPath)-9] = 0;
    string exec = execPath;
    cout << "Path to executables: " << exec << endl;
    
    if (exec == "")
        exec = "./";

    commandArg<string> iStringCmmd("-i","read fasta file");
    commandArg<string> iwormStringCmmd("-iworm","inchworm file", "");
    commandArg<string> oStringCmmd("-o","output directory");
    commandArg<bool> pairsStringCmmd("-paired", "paired-end reads are used.", false);
    commandArg<string> butterflyCmmd("-butterfly","butterfly executable", "../Butterfly/Butterfly.jar");
    commandArg<bool> skipCmmd("-skip","skip initial 2 steps", false);
    commandArg<bool> strandCmmd("-strand","strand-specific data", false);
    commandArg<bool> nobreakCmmd("-nobreak","skip breaking", false);
    commandArg<int> minCmmd("-min","minimum sequence length", 300);
    commandArg<int> cpuCmmd("-cpu","number of CPUs to use", 10);
    commandArg<int> distCmmd("-dist","size of a read pair insert", 350);
    commandArg<int> minDumpLenCmmd("-min_all","skip components for which all seqs < min_all", 110);
    commandArg<bool> buttCmmd("-run_butterfly","runs butterfly locally", false);
    commandArg<int>  kmerSizeCmmd("-kmer_size", "size of k-mer overlap length", KMER_SIZE);
    commandArg<int>  weldmerSizeCmmd("-weldmer_size", "size of weldmer to be used by GraphFromFasta", 48);
    commandArg<long> maxReadsCmd("-max_reads", "max number of reads to map to each graph", -1);
    commandArg<long> mmrCmmd("-max_mem_reads","Maximum number of reads to load into memory", -1);
    commandArg<int>  pctReadMap("-min_pct_read_mapping", "minimum percent of a read's kmers that must map to an inchworm bundle", 0);
    commandArg<int> minGlueCmmd("-min_glue", "minimum read support for glue in GraphToFasta", 2);
    commandArg<double> glueFactorCmmd("-glue_factor", "fraction of max (iworm pair coverage) for read glue support", 0.05);
    commandArg<double> minIsoRatioCmmd("-min_iso_ratio", "min ratio of (iworm pair coverage) for join", 0.05);
    commandArg<bool> debugCmmd("-debug", "verbosely describes operations", false);
    commandArg<string> sortBufferSizeCmmd("-sort_buffer_size","size of memory buffer to use in all sort calls", "2G");
    commandArg<bool> no_cleanupCmmd ("-no_cleanup", "retain input files on success", false);
    commandArg<string> readsForPairsCmmd("-reads_for_pairs", "reads fasta file to use for pairing analysis", "");
    commandArg<bool> noPairLinksCmmd("-no_pair_links", "ignore pair link info in clustering", false);
    commandArg<bool> bowtieCompCmmd("-bowtie_comp", "use Bowtie2 for readsToTranscripts mapping", false);
    commandArg<string> sortExecCmmd("-sort_exec", "sort command to use, ie. sort --parallel=${Nthreads}", sort_exec);
    commandArg<bool> debugAllWeldCmmd("-debug_weld_all", "cluster all contigs together, debugging only", false);

    
    commandLineParser P(argc,argv);
    P.SetDescription("Assemble transcriptomes from reads.");
    P.registerArg(iStringCmmd);
    P.registerArg(iwormStringCmmd);
    P.registerArg(butterflyCmmd);
    P.registerArg(oStringCmmd);
    P.registerArg(pairsStringCmmd);
    P.registerArg(skipCmmd);
    P.registerArg(strandCmmd);
    P.registerArg(nobreakCmmd);
    P.registerArg(minCmmd);
    P.registerArg(cpuCmmd);
    P.registerArg(buttCmmd);
    P.registerArg(distCmmd);
    P.registerArg(minDumpLenCmmd);
    P.registerArg(kmerSizeCmmd);
    P.registerArg(weldmerSizeCmmd);
    P.registerArg(maxReadsCmd);
    P.registerArg(mmrCmmd);
    P.registerArg(pctReadMap);
    P.registerArg(minGlueCmmd);
    P.registerArg(glueFactorCmmd);
    P.registerArg(minIsoRatioCmmd);
    P.registerArg(debugCmmd);
    P.registerArg(sortBufferSizeCmmd);
    P.registerArg(no_cleanupCmmd);
    P.registerArg(readsForPairsCmmd);
    P.registerArg(noPairLinksCmmd);
    P.registerArg(bowtieCompCmmd);
    P.registerArg(sortExecCmmd);
    P.registerArg(debugAllWeldCmmd);

    P.parse();

    cerr << "-------------------" << endl 
         << "---- Chrysalis ----" << endl
         << "-------------------" << endl << endl;
    

    string readString = P.GetStringValueFor(iStringCmmd);
    string readsForPairs = P.GetStringValueFor(readsForPairsCmmd);
    sort_exec = P.GetStringValueFor(sortExecCmmd);
    
    if (readsForPairs.size() == 0) {
        readsForPairs = readString;
    }
    
    string iwormString = P.GetStringValueFor(iwormStringCmmd);
    string outDir = P.GetStringValueFor(oStringCmmd);
    string butterflyExec = P.GetStringValueFor(butterflyCmmd);
    string sortBufferSizeString = P.GetStringValueFor(sortBufferSizeCmmd);
    

    bool PAIRED_READS_MODE = P.GetBoolValueFor(pairsStringCmmd);
    bool bSkip = P.GetBoolValueFor(skipCmmd);
    bool sStrand = P.GetBoolValueFor(strandCmmd);
    int minDumpLen = P.GetIntValueFor(minDumpLenCmmd);
    int minLen = P.GetIntValueFor(minCmmd);
    int nCPU = P.GetIntValueFor(cpuCmmd);
    int pairDist = P.GetIntValueFor(distCmmd);
    bool NO_PAIR_LINKS = P.GetBoolValueFor(noPairLinksCmmd);
    bool BOWTIE_COMP = P.GetBoolValueFor(bowtieCompCmmd);
     
    NO_CLEANUP = P.GetBoolValueFor(no_cleanupCmmd);
    __DEBUG_WELD_ALL = P.GetBoolValueFor(debugAllWeldCmmd);
    
    bool bBreak = true;
    if (P.GetBoolValueFor(nobreakCmmd))
        bBreak = false;
    
    double glue_factor = P.GetDoubleValueFor(glueFactorCmmd);
    bool bButt = P.GetBoolValueFor(buttCmmd);
    
    long max_reads = P.GetLongValueFor(maxReadsCmd);    
    long max_mem_reads = P.GetLongValueFor(mmrCmmd);
    
    int min_pct_read_mapping = P.GetIntValueFor(pctReadMap);
    int min_glue = P.GetIntValueFor(minGlueCmmd);
    double min_iso_ratio = P.GetDoubleValueFor(minIsoRatioCmmd);
    int weldmer_size = P.GetIntValueFor(weldmerSizeCmmd);
    KMER_SIZE = P.GetIntValueFor(kmerSizeCmmd);
    


    bool DEBUG = P.GetBoolValueFor(debugCmmd);
    

    if (minDumpLen > minLen) {
        minDumpLen = minLen;
    }
    
    string command;
    
    command = "mkdir ";
    command += outDir;
    command += " 2>/dev/null"; //ignore any error msg
    system(command.c_str());
    

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Scaffolding of Inchworm contigs:
    // -run bowtie to align reads to inchworm contigs, and determine scaffolding support based on pairs
    if (PAIRED_READS_MODE && (! NO_PAIR_LINKS)) {
        // Run bowtie to map reads to inchworm transcripts

        string scaffolds_finished_file = "iworm_scaffolds.txt.finished";
        if (Exists(scaffolds_finished_file) && Exists("iworm_scaffolds.txt")) {
            cerr << "Warning, Scaffolds file: iworm_scaffolds.txt previously generated and being reused." << endl;
        }
        else {
                    
            //////////////////////////////////
            // Running Bowtie Directly Now //
            /////////////////////////////////

            cerr << "###########################################################" << endl
                 << "## Running Bowtie to map fragment ends to Inchworm Contigs" << endl
                 << "## to use pairing info for Chrysalis' clustering " << endl
                 << "###############################################################" << endl << endl;
            
            stringstream cmdstr;
            cmdstr << "ln -sf " << iwormString << " target.fa ";
            cerr << "CMD: " << cmdstr.str() << endl;
            Execute(cmdstr.str().c_str());

            cmdstr.str(""); // clear it out.
            if (Exists("target.fa.finished")){
                 cerr << "Warning, bowtie-build output already exists and will be re-used: target.fa*" << endl;
            }else{
                 cmdstr << "bowtie-build -q target.fa target";
                 cerr << "CMD: " << cmdstr.str() << endl;
                 Execute(cmdstr.str().c_str());
                 string finished_file = "target.fa.finished";
                 string complete_cmd = "touch " + finished_file;
                 Execute(complete_cmd.c_str());
            }

            string bowtie_sam_file = "bowtie.nameSorted.bam";            
            string bowtie_checkpoint = bowtie_sam_file + ".finished";
            cmdstr.str(""); // clear it out.
            if (Exists(bowtie_checkpoint)){
                cerr << "Warning, bowtie output already exists and will be re-used: " << bowtie_sam_file << endl;
            }else{
                 cmdstr << "bash -c \" set -o pipefail; bowtie -a -m 20 --best --strata --threads " << nCPU 
                        << " --chunkmbs 512 -q -S "
                        << " -f target " << readsForPairs 
                        << " | samtools view -F4 -Sb - | samtools sort -no - - > " << bowtie_sam_file << " \" ";
                 
                 cerr << "CMD: " << cmdstr.str() << endl;
                 Execute(cmdstr.str().c_str());
                 string complete_cmd = "touch " + bowtie_checkpoint;
                 Execute(complete_cmd.c_str());
            }
                        
            // Generate inchworm pair links:
            cmdstr.str(""); // clear it out.
            
            cmdstr << exec << "/../util/support_scripts/scaffold_iworm_contigs.pl " << bowtie_sam_file << " " << iwormString << " > iworm_scaffolds.txt";
            cerr << "CMD: " << cmdstr.str() << endl;
            
            Execute(cmdstr.str().c_str());

            string complete_scaffolds_cmd = "touch " + scaffolds_finished_file;
            Execute(complete_scaffolds_cmd.c_str());
            
            // cleanup
            if (! NO_CLEANUP) {
                // purge the alignments
                string rm_cmd = "rm -rf iworm_bowtie/ ";
                Execute(rm_cmd.c_str());
            }
            
        }
    }
    

    /////////////////////////////////////////////////////////////////////////////////////
    // GraphFromFasta:
    // Pool sequences (clustering of Inchworm contigs, join using read-support)
    // Create output file: GraphFromFasta.out
    
    cerr << "-- running Chryalis: GraphFromFasta --" << endl;


    string graphFromFastaCompleteFile = outDir + "/GraphFromIwormFasta.finished";
    
    stringstream cmdstr;

    cmdstr << exec << "GraphFromFasta -i " << iwormString
           << " -r " << readString
           << " -min_contig_length " << minLen
           << " -min_glue " <<  min_glue
           << " -glue_factor " << glue_factor
           << " -min_iso_ratio " <<  min_iso_ratio
           << " -t " << nCPU 
           << " -k " << KMER_SIZE
           << " -kk " << weldmer_size;
    
    if (sStrand) {
        cmdstr << " -strand ";
    }
    
    // cmdstr << " -report_welds ";
    

    if (PAIRED_READS_MODE && (! NO_PAIR_LINKS)) {
        
        cmdstr << " -scaffolding iworm_scaffolds.txt ";

    }
    
    if (__DEBUG_WELD_ALL) {

        cmdstr << " -debug_weld_all ";
    }
    

    cmdstr << " > ";
    string components_file = outDir + "/GraphFromIwormFasta.out";
    cmdstr << components_file;
    
    command = cmdstr.str();
    
    if (Exists(components_file) && Exists(graphFromFastaCompleteFile)) {
        cerr << "File: " << components_file << " already exists. Using existing file assuming resume mode" << endl << endl;
    }
    else {
        cout << "Running: " << command << endl;
        Execute(command.c_str());

        // make completion-marking file
        string complete_file_cmd = "touch " + graphFromFastaCompleteFile;
        Execute(complete_file_cmd.c_str());
        
    }
    //
    ///////////////////////////////////////////////////////////////////////////////////////
    
    

    ///////////////////////////////////////////////////////////////////////////////////////
    // Read components.out, create chrysalis/bundled.fasta which represents clustered inchworm sequences
    // Create de Bruijn graphs based on clustered Inchworm seqs (partitioned chrysalis/RawComps.\d+/comp\d+.raw.graph files) 

    cerr << "-- writing inchworm bundled.fasta and computing & partitioning component graph files for parallel processing." << endl;


    string bundledName = outDir + "/bundled_iworm_contigs.fasta";
    string bundledFinished = bundledName + ".finished";

    if (Exists(bundledName) && Exists(bundledFinished)) {
        cerr << "Chrysalis' inchworm contig bundling already completed. Reusing existing bundles." << endl;
        sleep(2);
    }
    else {
        
        //---- Note: This section of code requires that you have the stacksize set to unlimited. --------//
        
        string iString = components_file;
        //string oString = outDir + "/graph.out";
        
        string tmpName = outDir + "/tmp.fasta";
        //FILE * p = fopen(oString.c_str(), "w");
        //fclose(p);
        
        FlatFileParser parser;  
        parser.Open(iString);  // read components.out file
        
        int component_no = 0;
        
        ComponentFileMgr mgr(outDir);
        
        vecDNAVector tmpSeq;
        vector<int> component_ids;
        
        vector<string_vec> bundled;
        //bundled.reserve(1000000);
        
        //DNAVector separator;
        //separator.SetFromBases("X");
        string separator = "X";
        
        FILE * pOut = NULL;
        
        while (parser.ParseLine()) {
            if (parser.GetItemCount() == 0)
                continue;
            
            if (parser.AsString(0)[0] == '#') {
                continue; // ignore comments
            }
            
            if (parser.AsString(0) == "COMPONENT") {
                component_no = parser.AsInt(1);
                int num_iworm_contigs = parser.AsInt(2);
                
                tmpSeq.resize(0); // hold the iworm seqs corresponding to a single component
                
                while (parser.ParseLine()) {
                    
                    if (parser.GetItemCount() == 0)
                        continue;
                    if (parser.AsString(0)[0] == '#') {
                        continue; // ignore comments
                    }         
                    
                    
                    if (parser.AsString(0) == "END") {
                        break;	  
                    }
                    
                    const char * ppp = parser.AsString(0).c_str();
                    
                    if (ppp[0] == '>') {
                        DNAVector tmp;
                        tmpSeq.push_back(tmp);
                        continue;
                    }
                    
                    DNAVector app;
                    app.SetFromBases(parser.AsString(0));
                    //cerr << "adding: " << parser.AsString(0) << endl;
                    tmpSeq[tmpSeq.isize()-1] += app; // adding the iworm sequence to that iworm entry
                    
                    //fprintf(p, "%s\n", parser.Line().c_str());
                }
                //fclose(p);
                
                if (tmpSeq.isize() == 1 && tmpSeq[0].isize() < minLen) {
                    //cerr << "-discarding entry, too short." << endl;
                    continue;
                }
                
                bool bGood = false;
                for (int x=0; x<tmpSeq.isize(); x++) {
                    if (tmpSeq[x].isize() > minDumpLen) {
                        bGood = true;
                        break;
                    }
                }
                
                if (bGood) {
                    
                    vector<string> iworm_bundle;
                    for (int x = 0; x < tmpSeq.isize(); x++) {
                        iworm_bundle.push_back(get_seq_string(tmpSeq[x]));    
                        
                    }
                    
                    bundled.push_back(iworm_bundle);
                    component_ids.push_back(component_no);
                }
                
                tmpSeq.resize(0);
                
            }
            
        }
        
        write_iworm_bundle(bundledName, bundled, component_ids, mgr);
     
        
        // make completion-marking file
        string complete_file_cmd = "touch " + bundledFinished;
        Execute(complete_file_cmd.c_str());
        
    }
    

    /////////////////////////////////////////////////////////////////////
    // Make fasta files for each component...
    
    cerr << "-- mapping reads to chrysalis components." << endl;
    
    string readcounts_file = outDir + "/rcts.out";
    string readsToTranscriptsCompleteFile = outDir + "/readsToComponents.finished";
    string readsToTranscriptsOutputFile = outDir + "/readsToComponents.out";
    
    if (! Exists(readsToTranscriptsCompleteFile) ) {

        stringstream command;
        

        if( BOWTIE_COMP ){
          command << exec << "ReadsToComponents.pl ";
        } else {
          command << exec << "ReadsToTranscripts ";
        }
        command << " -i " << readString
                << " -f  " << bundledName
                << " -t " << nCPU
                << " -o " << outDir;
        if (sStrand) {
          command << " -strand ";
        }
        if (min_pct_read_mapping > 0) {
          command << " -p " << min_pct_read_mapping;
        }
        if (max_mem_reads > 0) {
          command << " -max_mem_reads " << max_mem_reads;
        }
        cout << "Running: " << command.str() << endl;
        Execute(command.str().c_str());

        struct stat filestatus;
        stat(readsToTranscriptsOutputFile.c_str(), &filestatus );
        cout << filestatus.st_size << " bytes\n";
        if (filestatus.st_size == 0) {
            cerr << "ERROR, " << readsToTranscriptsOutputFile << " has zero file size. Read mappings to inchworm components must have failed" << endl;
            exit(3);
        }

        string complete_cmd = "touch " + readsToTranscriptsCompleteFile;
        Execute(complete_cmd.c_str());
    }
    else {
        cerr << "File: " << readsToTranscriptsCompleteFile << " already exists. Using existing read/transcript mappings assuming resume mode." << endl << endl;
    }
    
    

    /*
    // readcounts_file created by the above.
    FlatFileParser readCount;  
    readCount.Open(readcounts_file);
    readCount.ParseLine();
    string numReads = readCount.Line();
    */  // dont need that anymore.
   
    
    // sort output by component identifier.
    string sortedReadsToTranscriptsOutputFile = readsToTranscriptsOutputFile + ".sort";
    string sortedReadsToTranscriptsFinishedFile = sortedReadsToTranscriptsOutputFile + ".finished";
    
    if (! Exists(sortedReadsToTranscriptsFinishedFile)) {

        cmdstr.str(""); // clear 
        cmdstr << sort_exec << " -T . -S " << sortBufferSizeString 
               << " -k 1,1n " << readsToTranscriptsOutputFile  << " > " << sortedReadsToTranscriptsOutputFile;
        
        cerr << "CMD: " << cmdstr.str() << endl;
        Execute(cmdstr.str().c_str());


        struct stat filestatus;
        stat(sortedReadsToTranscriptsOutputFile.c_str(), &filestatus );
        cout << filestatus.st_size << " bytes\n";
        if (filestatus.st_size == 0) {
            cerr << "ERROR, " << sortedReadsToTranscriptsOutputFile << " has zero file size. Sorting of read mappings to inchworm components must have failed" << endl;
            exit(3);
        }
        
        
        string complete_cmd = "touch " + sortedReadsToTranscriptsFinishedFile;
        Execute(complete_cmd.c_str());
    
        if (! NO_CLEANUP) {
            // now that we have the sorted version, no longer need the unsorted one.
            unlink(readsToTranscriptsOutputFile.c_str());
            
        }
        
    }
    
    //////////////////////////////////////////////////////////////////////
    // QuantifyGraph
    //////////////////////////////////////////////////////////////////////


    /*

    cerr << "-- writing quantifyGraph and butterfly commands for parallel processing." << endl;
    
    //svec<string> targets;
    string butterflyCommands = outDir + "/butterfly_commands";
    FILE * pButterfly = fopen(butterflyCommands.c_str(), "w");

    string quantifyGraphCommands = outDir + "/quantifyGraph_commands";
    FILE * pQGraphCmds = fopen(quantifyGraphCommands.c_str(), "w");

    ofstream component_listing_ofh;
    string component_listing_filename = outDir + "/component_file_listing.txt";
    component_listing_ofh.open(component_listing_filename.c_str());
    
    map<int,bool> pursue_component;

    for (int i=0; i<component_ids.size(); i++) {
        int component_no = component_ids[i];
        string graph = mgr.GetFileName(component_no, ".raw.graph", false);
        string fasta = mgr.GetFileName(component_no, ".raw.fasta", false);
        string finalgraph = mgr.GetFileName(component_no, ".out", false);
        //string fnialreads = mgr.GetFileName(k, ".finalgraph.reads");
        
        string component_file_basename = mgr.GetFileName(component_no, "", false);
        
        FILE * pFasta = fopen(fasta.c_str(), "r");
        if (pFasta == NULL)
            continue;
       
        fclose(pFasta);
        
        pursue_component[component_no] = true;
        
        char cwd[FILENAME_MAX];
        getcwd(cwd, sizeof(cwd) / sizeof(char));

        stringstream bfly_cmd;

        bfly_cmd << "java -jar " << butterflyExec << " -N " << numReads << " -L " << minLen << " -F " << pairDist;
        //bfly_cmd << " -C " << cwd << "/" << mgr.GetFileName(component_no, "", false);
        bfly_cmd << " -C " << mgr.GetFileName(component_no, "", false); // filenames now provided as fully qualified.
        
        if (NO_CLEANUP) {
            bfly_cmd << " --no_cleanup ";
        }


        fprintf(pButterfly, "%s\n", bfly_cmd.str().c_str());
        

        stringstream qgraph_cmd;
        //qgraph_cmd << exec << "QuantifyGraph -g " << cwd << "/" << graph 
        //           << " -o " << cwd << "/" << finalgraph << " -i " << cwd << "/" << fasta;
        qgraph_cmd << exec << "QuantifyGraph -g " << graph 
                   << " -o " << finalgraph << " -i " << fasta;
        
        if (sStrand) {
            qgraph_cmd << " -strand";
        }
        

        if (max_reads > 0) {
            char max_reads_string[256];
            sprintf(max_reads_string, " -max_reads %li ", max_reads);
            // command += max_reads_string;
            qgraph_cmd << " -max_reads " << max_reads;
        }
        
        if (NO_CLEANUP) {
            qgraph_cmd << " -no_cleanup ";
        }
        
        
        fprintf(pQGraphCmds, "%s\n", qgraph_cmd.str().c_str());

        component_listing_ofh << component_file_basename << endl;




    }
    fclose(pButterfly);
    fclose(pQGraphCmds);
    component_listing_ofh.close();


    */

    // write the deBruijn graphs:
    //write_deBruijn_graphs(bundled, component_ids, mgr, outDir, exec, sStrand, pursue_component);

    write_deBruijn_graphs(bundledName, exec, sStrand, nCPU);
    
    
    return(0);

}
