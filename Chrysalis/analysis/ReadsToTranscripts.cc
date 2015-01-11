/* -*- mode: C++; c-basic-offset: 4; indent-tabs-mode: nil -*- */
#include <string>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() (1)
#define omp_get_num_threads() (1)
#define omp_get_thread_num()  (0)
#endif

#include "base/CommandLineParser.h"
#include "analysis/DNAVector.h"
#include "analysis/NonRedKmerTable.h"

//#include "analysis/CompMgr.h"
#include <queue>
#include <map>
#include <algorithm>


#define MAX_OPEN_FILES 1000

extern "C"
{
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
}

class Assignment
{
public:
    Assignment() {
        m_comp = -1;
        m_read = -1;
    }

    void SetComp(int c) {
        m_comp = c;
    }
    void SetRead(int r) {
        m_read = r;
    }

    int Read() const {return m_read;}
    int Comp() const {return m_comp;}

    bool operator < (const Assignment & a) const {
        return (m_comp < a.m_comp);
    }

private:
    int m_comp;
    int m_read;
};



int main(int argc,char** argv)
{


    commandArg<string> aStringCmmd("-i","reads fasta");
    commandArg<string> fStringCmmd("-f","fasta input file (concatenated flat components)");
    commandArg<string> bStringCmmd("-o","output file");
    commandArg<long> mmrCmmd("-max_mem_reads","Maximum number of reads to load into memory", -1);
    commandArg<bool> strandCmmd("-strand","strand specific data", false);
    commandArg<int> pctReadMapCmmd("-p", "percent of read kmers require mapping to component", 0);
    commandArg<bool> verboseCmmd("-verbose", "prints more status info", false);
    commandArg<int> numThreadsCmmd("-t", "number of threads (default: env OMP_NUM_THREADS)", 0);

    commandLineParser P(argc,argv);
    P.SetDescription("Assigns reads to graph components.");
    P.registerArg(aStringCmmd);
    P.registerArg(fStringCmmd);
    P.registerArg(bStringCmmd);
    P.registerArg(mmrCmmd);
    P.registerArg(strandCmmd);
    P.registerArg(pctReadMapCmmd);
    P.registerArg(verboseCmmd);
    P.registerArg(numThreadsCmmd);

    P.parse();

    cerr << "-------------------------------------------" << endl
         << "---- Chrysalis: ReadsToTranscripts --------" << endl
         << "-- (Place reads on Inchworm Bundles) ------" << endl
         << "-------------------------------------------" << endl << endl;
    

    string aString = P.GetStringValueFor(aStringCmmd); // reads.fa
    string readsToTranscriptsOutputFile = P.GetStringValueFor(bStringCmmd); // out directory
    string fString = P.GetStringValueFor(fStringCmmd); // inchworm bundled.fasta
    long max_mem_reads = P.GetLongValueFor(mmrCmmd);
    bool bStrand = P.GetBoolValueFor(strandCmmd);  // strand-specific data
    int pctReadMapRequired = P.GetIntValueFor(pctReadMapCmmd); 
    int num_threads = P.GetIntValueFor(numThreadsCmmd);
    
    bool VERBOSE = P.GetBoolValueFor(verboseCmmd);
    
    vecDNAVector dna;

    if(max_mem_reads > 0){
        cout << "Setting maximum number of reads to load in memory to " << max_mem_reads << endl;
    } else {
        max_mem_reads = 2147483647; // max int
    }


    if (num_threads > 0) {
        cerr << "-setting num threads to: " << num_threads << endl;
        omp_set_num_threads(num_threads);
    }

    cerr << "Reading bundled inchworm contigs... " << endl;
    //ML: use 1M as read buffer size
    dna.Read(fString, false, false, true, 1000000);
    cerr << "done!" << endl;
    //test.ReverseComplement();

    int k = 25;

    //ComponentFileMgr compMgr(bString);
    multimap<int,int> componentReadMap;

    unsigned long readCount = 0;
    string readCountFile = readsToTranscriptsOutputFile + ".rcts.out";


    //-------------------------------------------------//
    //---- Build kmer index for Iworm Bundles ---------//
    //-------------------------------------------------//
    

    NonRedKmerTable kt(k);
    kt.SetUp(dna, true);
    kt.SetAllCounts(-1);

    map<int,int> component_number_mapping;

    cerr << "Assigning kmers to Iworm bundles ... ";
    #pragma omp parallel for schedule(guided,100)
    for (int i=0; i< dna.isize(); i++) {
        const DNAVector & d = dna[i];
        string bundle_acc = dna.Name(i);
        string component_no_string = bundle_acc.substr(3); // remove >s_ prefix
        int component_no = atoi(component_no_string.c_str());
        
        #pragma omp critical
        {
            component_number_mapping[i] = component_no;
            
            //cerr << "\rBundle name: [" << i << "] = " << bundle_acc << ", so component_no = " << component_no << "    ";
        }


        for (int j=0; j<=d.isize()-k; j++) {
            kt.SetCount(d, j, i);  // not really counting kmers here, but instead labeling kmers according to the iworm bundles they correspond to.
        }
    }
    cerr << "done!" << endl;


    // ------------------------------------------------//
    // ------- Assign reads to Iworm Bundles ----------//
    // ------------------------------------------------//

    //ML: this implementation does not make use of class DNAVector for the reads

    map<int,bool> seen_component; // open files on first write, append on later writes.

    vector<int> read_pct_mapping_info;
    read_pct_mapping_info.resize(100000, 0);

    DNAStringStreamFast fastaReader;
    fastaReader.ReadStream(aString);


    //string readsToTranscriptsOutputFile = bString + "/readsToComponents.out"; // really mapping reads to components rather than transcripts, using spectral alignment
    int outFile = open(readsToTranscriptsOutputFile.c_str(), O_WRONLY | O_CREAT | O_TRUNC , 0666);
        
    cerr << "Processing reads:" << endl;
    unsigned long total_reads_read = 0;
    int reads_read = 0;
    do {
        vector<string> readSeqVector;
        vector<string> readNameVector;
        
        cerr << " reading another " << max_mem_reads << "... ";
        for(int i=0;i<max_mem_reads;i++)
        {
            if(!fastaReader.NextToVector(readSeqVector, readNameVector)) break;
        }

        reads_read = readSeqVector.size();
        if (reads_read > 0) {
            cerr << "done.  Read " << reads_read << " reads." << endl;
        } else {
            // reached EOF
            cerr << "finished reading reads" << endl;
        }
        
        componentReadMap.clear();
    
        #pragma omp parallel for schedule(guided,100)
        for (int i=0; i < (int)readSeqVector.size(); i++) {
            
            string d = readSeqVector[i];
            // ensure upper case
            transform(d.begin(), d.end(), d.begin(), ::toupper);
            
            // map kmers of read to Iworm bundles
            svec<int> comp;
            comp.reserve(4000);
            int num_kmer_pos = d.size()-k + 1;
            for (int j=0; j<=(int)d.size()-k; j++) {
                int c = kt.GetCountReal(d, j); // the iworm bundle containing read kmer
                if (c >= 0)
                    comp.push_back(c);
            }
            if (!bStrand) {
                string dd = d;
                DNAVector::ReverseComplement(dd);
                for (int j=0; j<=(int)dd.size()-k; j++) {
                    int c = kt.GetCountReal(dd, j);
                    if (c >= 0)
                        comp.push_back(c);
                }
            }

            // find the iworm bundle with the most kmer hits
            Sort(comp);
            int best = -1;
            int max = 0;
            int run = 0;
            for (int j=1; j<comp.isize(); j++) {
                if (comp[j] != comp[j-1] || j+1 == comp.isize()) {
                    if (run > max) {
                        max = run;
                        best = comp[j-1];
                    }
                    run = 0;
                } else {
                    run++;
                }
            }
            
            //cout << "Read " << i << " maps to " << best << " with " << max << " k-mers" << endl;
            int pct_read_mapped = (int) ((float)max/num_kmer_pos*100 + 0.5);
            //ML: the following calculates the same, but without floating point precision issues:
            //int pct_read_mapped = ( 100*max + num_kmer_pos/2 ) / num_kmer_pos;
            if(best != -1 && pct_read_mapped >= pctReadMapRequired) {
                // Note: multiple threads shouldn't try to reorder the tree
                // at the same time (which can happen during an insert)
                #pragma omp critical
                {
                    readCount++;
                    componentReadMap.insert(pair<int,int>(best,i));  // i = read_index, best = best_iworm_bundle
                    
                    // store pct read mapping info for later reporting
                    if (i > (int) read_pct_mapping_info.size()-1) {
                        read_pct_mapping_info.resize(i+100000);
                        // cerr << "resizing vec to " << i << " + 100000 " << endl;
                    }
                    read_pct_mapping_info[i] = pct_read_mapped;
                    
                }
            }
            else {
                // no mapping for read
                if (VERBOSE)
                    cerr << "WARNING: No component mapping for read: " << readNameVector[i] << " : " << readSeqVector[i] << endl;
            }
           
        } // end of read assignment to components for this round of streaming
        
        
        if (reads_read > 0) {
            total_reads_read += reads_read;
            cerr << "[" << total_reads_read << "] reads analyzed for mapping." << endl;
        }

        //-------------------------------------------------------
        // Write reads to files based on components mapped to.

        //cerr << "\nWriting to files... ";
        // convert sorted map to vectors for threaded processing
        vector<int> mapComponents; // list of component_ids
        vector<vector<int> > componentReads; // list of read-vectors that correspond to the component_ids above (synched by index)
        multimap<int,int>::iterator it;
        int lastComponent = -1;
        vector<int> tmpReadIDs; // temporary hold for the reads corresponding to the last component
        for(it = componentReadMap.begin(); it != componentReadMap.end(); it++){
            if(it->first != lastComponent){  // start new component, first store previous component info.
                if(tmpReadIDs.size() > 0){
                    mapComponents.push_back(lastComponent);
                    componentReads.push_back(tmpReadIDs);
                }
                lastComponent = it->first;
                tmpReadIDs.clear();
            }
            tmpReadIDs.push_back(it->second);
        }
        if(tmpReadIDs.size() > 0){ // get the last component info stored.
            mapComponents.push_back(lastComponent);
            componentReads.push_back(tmpReadIDs);
        }


        // write out to files in map order
         
        #pragma omp parallel for schedule(static,1)
        for(unsigned int i = 0; i < mapComponents.size(); i++){
            
            // each thread accesses a unique component number, so no competition between threads in writing to component-based files.
            
            int iworm_bundle_index = mapComponents[i];
            int iworm_component_no = component_number_mapping[iworm_bundle_index];

            /*  changed to write to a single file - bhaas  July 3, 2013
            
            string name = compMgr.GetFileName(iworm_component_no, ".raw.fasta", false);
            int outFile;
            
            //ML: use systm calls open/write/close, buffer explicitly (faster than C++ buffered I/O)
            if(seen_component.find(iworm_component_no) == seen_component.end()) {
                // first write
                if (VERBOSE)
                    cerr << "First write to component: " << iworm_component_no << endl;
                outFile = open(name.c_str(), O_WRONLY | O_CREAT | O_TRUNC , 0666);
                #pragma omp critical
                seen_component[iworm_component_no] = true;
            } 
            else {
                // subsequent write as append
                if (VERBOSE)
                    cerr << "Second write as append to component: " << iworm_component_no << endl;
                outFile = open(name.c_str(), O_WRONLY | O_APPEND );
            }
            if(outFile<0) cerr << "cannot open " << name << endl;
            
            */


            //ML: put everthing to write into buf
            string tmpName;
            char component_id_string[8];
            sprintf(component_id_string, "%d", iworm_component_no);
            string buf;
            for(unsigned int j = 0; j < componentReads[i].size(); j++){
                buf += component_id_string;
                buf += "\t";
                int read_index = componentReads[i][j];
                DNAStringStreamFast::formatReadNameString(readNameVector[read_index], tmpName);
                buf += tmpName;
                char foo[8];
                sprintf(foo, "%d%%\t", read_pct_mapping_info[read_index]);
                buf += "\t";
                buf += foo;
                buf += readSeqVector[read_index];
                buf += "\n";
            }
            
            // write mapping to file, avoid thread collisions
            #pragma omp critical
            if(write(outFile,buf.c_str(),buf.size())<0)
            {
                cerr <<  "error writing file "  << readsToTranscriptsOutputFile  << ": " <<  strerror(errno)  << endl;
                exit(0);
            }
            //close(outFile);
        }

        if(mapComponents.size()>0)
          cerr << "[" << mapComponents.size() << "] components written." << endl;
        
        //cerr << "done\n";
        //clear out read component mappings
    } while (reads_read > 0);
    
    close(outFile);
    
    cerr << "Done" << endl;

    FILE * pReadCount = fopen(readCountFile.c_str(), "w");
    fprintf(pReadCount, "%lu\n", readCount/2);
    fclose(pReadCount);

    return 0;
}

