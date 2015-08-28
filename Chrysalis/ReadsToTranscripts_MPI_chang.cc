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

#include "analysis/CompMgr.h"
#include <queue>
#include <map>

#ifdef MPI_ENABLED
#include<mpi.h>
#endif

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

#include<sys/time.h>
static struct timeval start,end;
void timer_start(){

  gettimeofday(&start,NULL);
}


double timer_stop(){
  gettimeofday(&end,NULL);
  double time_taken = ((end.tv_usec-start.tv_usec) + 1000000*(end.tv_sec - start.tv_sec));
 time_taken = time_taken/1000;
  return time_taken;
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

  int rank,numranks;
#ifdef MPI_ENABLED
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
  rank=0;
  numranks=1;
#endif

    commandArg<string> aStringCmmd("-i","reads fasta");
    commandArg<string> fStringCmmd("-f","fasta input file (concatenated flat components)");
    commandArg<string> bStringCmmd("-o","output directory");
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
    string bString = P.GetStringValueFor(bStringCmmd); // out directory
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

    ComponentFileMgr compMgr(bString);
    multimap<int,int> componentReadMap;

    unsigned long readCount = 0;
    string readCountFile = bString + "/rcts.out";


    //-------------------------------------------------//
    //---- Build kmer index for Iworm Bundles ---------//
    //-------------------------------------------------//
    

    NonRedKmerTable kt(k);
    kt.SetUp(dna, true);
    kt.SetAllCounts(-1);

    map<int,int> component_number_mapping;

    cerr << "Assigning kmers to Iworm bundles ... on rank " << rank << endl;
    #pragma omp parallel for schedule(guided,100)
    for (size_t i=0; i< dna.size(); i++) {
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

#ifdef MPI_ENABLED
  if(rank==0) {
    fastaReader.ReadStream(aString);
  }
#else
    fastaReader.ReadStream(aString);
#endif

    string readsToTranscriptsOutputFile = bString + "/readsToComponents.out"; // really mapping reads to components rather than transcripts, using spectral alignment

int outFile;
#ifdef MPI_ENABLED
   if(rank==0) {
     outFile = open(readsToTranscriptsOutputFile.c_str(), O_WRONLY | O_CREAT | O_TRUNC , 0666);
   }
#else
     outFile = open(readsToTranscriptsOutputFile.c_str(), O_WRONLY | O_CREAT | O_TRUNC , 0666);
#endif 

  cerr << "Processing reads:" << endl;
  unsigned long total_reads_read = 0;
  int reads_read;
  do {

#ifdef MPI_ENABLED

  reads_read = 0;

  int numDNAs;
  int numName;
  int numVector;
  char* DNAs_chars;
  char* Name_chars;
  int*  DNAs_size;
  int*  Name_size;

  if(rank==0) {

	vector<string> readSeqVector;
        vector<string> readNameVector;

	cerr << " reading another " << max_mem_reads << "... " << endl;		
	for(int crank=1; crank<numranks; crank++) {

	  readSeqVector.clear();
          readNameVector.clear();

	  for(int i=0;i<max_mem_reads;i++)
          {
            if(!fastaReader.NextToVector(readSeqVector, readNameVector)) break;
          }
	
	  if( readSeqVector.size() > 0 ) {

	  	reads_read += (int)readSeqVector.size();

		DNAs_size = (int*)malloc(sizeof(int)*max_mem_reads);
        	Name_size = (int*)malloc(sizeof(int)*max_mem_reads);

		numDNAs = 0;
		numName = 0;
		for(int j=0; j<(int)readSeqVector.size(); j++) {
			const string d = readSeqVector[j];
			numDNAs += (int)d.size();
			const string Nm = readNameVector[j];
			numName += (int)Nm.size();
			DNAs_size[j] = (int)d.size(); 
			Name_size[j] = (int)Nm.size();	
		}

		numVector  = (int)readSeqVector.size();
		DNAs_chars = (char*)malloc(sizeof(char)*numDNAs);
		Name_chars = (char*)malloc(sizeof(char)*numName);

		int ld = 0;
		int ln = 0;
		for(int j=0; j<(int)readSeqVector.size(); j++) {
			const string d = readSeqVector[j];
			for(int k=0; k<(int)d.size(); k++) {
				DNAs_chars[ld] = d[k]; ld++;
			}
			const string Nm = readNameVector[j];
			for(int k=0; k<(int)Nm.size(); k++) { 
				Name_chars[ln] = Nm[k];	ln++;
			}
		}
			
		MPI_Send(&numVector, 1, MPI_INT, crank, 1, MPI_COMM_WORLD);	
		MPI_Send(&numDNAs,   1, MPI_INT, crank, 1, MPI_COMM_WORLD);
		MPI_Send(&numName,   1, MPI_INT, crank, 1, MPI_COMM_WORLD);

		MPI_Send(DNAs_chars, numDNAs,   MPI_CHAR, crank, 1, MPI_COMM_WORLD);
		MPI_Send(Name_chars, numName,   MPI_CHAR, crank, 1, MPI_COMM_WORLD);
		MPI_Send(DNAs_size,  numVector, MPI_INT,  crank, 1, MPI_COMM_WORLD);
		MPI_Send(Name_size,  numVector, MPI_INT,  crank, 1, MPI_COMM_WORLD);

		free(DNAs_chars);
		free(Name_chars);
		free(DNAs_size);
		free(Name_size);
	
	  } else {
		numVector = 0;
		MPI_Send(&numVector, 1, MPI_INT, crank, 1, MPI_COMM_WORLD);
	  }

        }

        if (reads_read > 0) {
            cerr << "done" << endl;
        } else {
            // reached EOF
            cerr << "finished reading reads "  << endl;
	  for(int crank=1; crank<numranks; crank++) {
	    numVector  = 0;
	    MPI_Send(&numVector, 1, MPI_INT, crank, 1, MPI_COMM_WORLD);
	  }
        }
  }
#else
        vector<string> readSeqVector;
        vector<string> readNameVector;
        
        cerr << " reading another " << max_mem_reads << "... ";
        for(int i=0;i<max_mem_reads;i++)
        {
            if(!fastaReader.NextToVector(readSeqVector, readNameVector)) break;
        }

        reads_read = readSeqVector.size();
        if (reads_read > 0) {
            cerr << "done" << endl;
        } else {
            // reached EOF
            cerr << "finished reading reads" << endl;
        }
#endif
       
#ifdef MPI_ENABLED

  int* read_pct_mapping_send;
  int* componentMap_send_first;
  int* componentMap_send_second;
  int rank_readCount;
  int numRPMI;

if(rank>0) {

    read_pct_mapping_info.resize(100000, 0);

    MPI_Recv(&numVector, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	

    if(numVector>0) {

	MPI_Recv(&numDNAs,   1,        MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&numName,   1,        MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	DNAs_chars = (char*)malloc(sizeof(char)*numDNAs);
        Name_chars = (char*)malloc(sizeof(char)*numName);

	MPI_Recv(DNAs_chars, numDNAs,  MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(Name_chars, numName,  MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        DNAs_size = (int*)malloc(sizeof(int)*max_mem_reads);
        Name_size = (int*)malloc(sizeof(int)*max_mem_reads);

        MPI_Recv(DNAs_size, numVector, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(Name_size, numVector, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	vector<string> readSeqVector;
        vector<string> readNameVector; 

	int td=0;
	int tn=0;
	for(int i=0; i<numVector; i++) {
		const string s(&DNAs_chars[td], DNAs_size[i]);	
		readSeqVector.push_back(s);

		const string ss(&Name_chars[tn], Name_size[i]);
		readNameVector.push_back(ss);

	        td += DNAs_size[i];
		tn += Name_size[i];
	}
	
	rank_readCount = 0;	
	componentReadMap.clear();	

	#pragma omp parallel for schedule(guided,100)
        for (int i=0; i < (int)readSeqVector.size(); i++) {

            const string d = readSeqVector[i];

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


	//  cout << "Read " << i << " maps to " << best << " with " << max << " k-mers" << endl;
	    int pct_read_mapped = (int) ((float)max/num_kmer_pos*100 + 0.5);
 
	    //ML: the following calculates the same, but without floating point precision issues:
	    //int pct_read_mapped = ( 100*max + num_kmer_pos/2 ) / num_kmer_pos;
	    if(best != -1 && pct_read_mapped >= pctReadMapRequired) {
		// Note: multiple threads shouldn't try to reorder the tree
		// at the same time (which can happen during an insert)
		#pragma omp critical
                {
                    rank_readCount++;
                    componentReadMap.insert(pair<int,int>(best,i));  // i = read_index, best = best_iworm_bundle

		    // store pct read mapping info for later reporting
		    if (i > (int) read_pct_mapping_info.size()-1) {
                        read_pct_mapping_info.resize(i+100000);
			// cerr << "resizing vec to " << i << " + 100000 " << endl;
		    }
		    read_pct_mapping_info[i] = pct_read_mapped;

                }
            }

        }  // end of read assignment to components for this round of streaming

	numRPMI = read_pct_mapping_info.size();
	
	read_pct_mapping_send    = (int*)malloc(sizeof(int)*rank_readCount);
	componentMap_send_first  = (int*)malloc(sizeof(int)*rank_readCount);
	componentMap_send_second = (int*)malloc(sizeof(int)*rank_readCount);	

	int tnum = 0;
	multimap<int,int>::iterator it;
	for(it = componentReadMap.begin(); it != componentReadMap.end(); it++){
		componentMap_send_first[tnum] = it->first;
		componentMap_send_second[tnum] = it->second;
		read_pct_mapping_send[tnum] = read_pct_mapping_info[tnum];
		tnum++;
	}

	MPI_Send(&rank_readCount,  1,  MPI_INT,  0, 2, MPI_COMM_WORLD);

	MPI_Send(&numVector, 1,        MPI_INT,  0, 2, MPI_COMM_WORLD);
	MPI_Send(&numDNAs,   1,        MPI_INT,  0, 2, MPI_COMM_WORLD);
        MPI_Send(&numName,   1,        MPI_INT,  0, 2, MPI_COMM_WORLD);
	MPI_Send(&numRPMI,   1,        MPI_INT,  0, 2, MPI_COMM_WORLD);
        MPI_Send(DNAs_chars, numDNAs,  MPI_CHAR, 0, 2, MPI_COMM_WORLD);
        MPI_Send(Name_chars, numName,  MPI_CHAR, 0, 2, MPI_COMM_WORLD);
        MPI_Send(DNAs_size, numVector, MPI_INT,  0, 2, MPI_COMM_WORLD);
        MPI_Send(Name_size, numVector, MPI_INT,  0, 2, MPI_COMM_WORLD);

	MPI_Send(componentMap_send_first,  rank_readCount, MPI_INT,  0, 2, MPI_COMM_WORLD);
	MPI_Send(componentMap_send_second, rank_readCount, MPI_INT,  0, 2, MPI_COMM_WORLD);
	MPI_Send(read_pct_mapping_send,    rank_readCount, MPI_INT,  0, 2, MPI_COMM_WORLD);

	free(DNAs_size);
        free(Name_size);
        free(DNAs_chars);
        free(Name_chars);

	free(read_pct_mapping_send);
	free(componentMap_send_first);
	free(componentMap_send_second);

    } else {
	rank_readCount = 0;
	MPI_Send(&rank_readCount,  1, MPI_INT,  0, 2, MPI_COMM_WORLD );
    }	

}

#else
        componentReadMap.clear();
    
        #pragma omp parallel for schedule(guided,100)
        for (int i=0; i < (int)readSeqVector.size(); i++) {
            
            const string d = readSeqVector[i];

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
           
        } // end of read assignment to components for this round of streaming
        
#endif 

#ifdef MPI_ENABLED
if(rank==0) {
        if (reads_read > 0) {
            total_reads_read += reads_read;
            cerr << "[" << total_reads_read << "] reads analyzed for mapping." << endl;
        }
}
#else
	if (reads_read > 0) {
            total_reads_read += reads_read;
            cerr << "[" << total_reads_read << "] reads analyzed for mapping." << endl;
        }
#endif

#ifdef MPI_ENABLED
	//-------------------------------------------------------
	// Write reads to files based on components mapped to.

if(rank==0) {

  for(int crank=1; crank<numranks; crank++) {

    MPI_Recv(&rank_readCount, 1, MPI_INT,  crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if( rank_readCount > 0) {

	readCount += rank_readCount;

	MPI_Recv(&numVector, 1,        MPI_INT, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&numDNAs,   1,        MPI_INT, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&numName,   1,        MPI_INT, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&numRPMI,   1,        MPI_INT, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        DNAs_chars = (char*)malloc(sizeof(char)*numDNAs);
        Name_chars = (char*)malloc(sizeof(char)*numName);

        MPI_Recv(DNAs_chars, numDNAs,  MPI_CHAR, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(Name_chars, numName,  MPI_CHAR, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        DNAs_size = (int*)malloc(sizeof(int)*max_mem_reads);
        Name_size = (int*)malloc(sizeof(int)*max_mem_reads);

        MPI_Recv(DNAs_size, numVector, MPI_INT, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(Name_size, numVector, MPI_INT, crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	vector<string> readSeqVector;
        vector<string> readNameVector;

        int td=0;
        int tn=0;
        for(int i=0; i<numVector; i++) {
                const string s(&DNAs_chars[td], DNAs_size[i]);
                readSeqVector.push_back(s);

                const string ss(&Name_chars[tn], Name_size[i]);
                readNameVector.push_back(ss);

                td += DNAs_size[i];
                tn += Name_size[i];
        }

	free(DNAs_size);
        free(Name_size);
        free(DNAs_chars);
        free(Name_chars);

	read_pct_mapping_send    = (int*)malloc(sizeof(int)*rank_readCount);
        componentMap_send_first  = (int*)malloc(sizeof(int)*rank_readCount);
        componentMap_send_second = (int*)malloc(sizeof(int)*rank_readCount);

        MPI_Recv(componentMap_send_first,  rank_readCount, MPI_INT,  crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(componentMap_send_second, rank_readCount, MPI_INT,  crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(read_pct_mapping_send,    rank_readCount, MPI_INT,  crank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	componentReadMap.clear();	
	read_pct_mapping_info.clear();
	read_pct_mapping_info.resize(numRPMI);

	for(int count=0; count<rank_readCount; count++) {
		componentReadMap.insert(pair<int,int>(componentMap_send_first[count],componentMap_send_second[count]));
		read_pct_mapping_info[ componentMap_send_second[count] ] = read_pct_mapping_send[count];
	}

	free(read_pct_mapping_send);
        free(componentMap_send_first);
        free(componentMap_send_second);

	vector<int> mapComponents;
	vector<vector<int> > componentReads;
	multimap<int,int>::iterator it;
	int lastComponent = -1;
	vector<int> tmpReadIDs;
	for(it = componentReadMap.begin(); it != componentReadMap.end(); it++){	
	    if(it->first != lastComponent){
		if(tmpReadIDs.size() > 0){
                    mapComponents.push_back(lastComponent);
                    componentReads.push_back(tmpReadIDs);
                }
                lastComponent = it->first;
                tmpReadIDs.clear();
            }
            tmpReadIDs.push_back(it->second);
        }
	if(tmpReadIDs.size() > 0){
	    mapComponents.push_back(lastComponent);
            componentReads.push_back(tmpReadIDs);
        }

	#pragma omp parallel for schedule(static,1)
        for(unsigned int i = 0; i < mapComponents.size(); i++){

	    int iworm_bundle_index = mapComponents[i];
            int iworm_component_no = component_number_mapping[iworm_bundle_index];

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

	    #pragma omp critical
            if(write(outFile,buf.c_str(),buf.size())<0)
            {
                cerr <<  "error writing file "  << readsToTranscriptsOutputFile << "on rank " << rank  << ": " <<  strerror(errno)  << endl;
                exit(0);
            }

	}

	if(mapComponents.size()>0)
          cerr << "[" << mapComponents.size() << "] components written. on rank " << crank << endl;	
    } 

  }

}

#else
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

#endif 

#ifdef MPI_ENABLED
if( rank == 0 ) {
  for(int crank=1; crank<numranks; crank++)
    MPI_Send(&reads_read, 1, MPI_INT,  crank, 3, MPI_COMM_WORLD);
} else {
    MPI_Recv(&reads_read, 1, MPI_INT,  0,     3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
#endif
        //cerr << "done\n";
        //clear out read component mappings
    } while (reads_read > 0);
  
#ifdef MPI_ENABLED
if(rank==0) {
    close(outFile);
    cerr << "Done" << endl;
}
#else
    close(outFile);
    cerr << "Done" << endl;
#endif

#ifdef MPI_ENABLED
if(rank==0) {
    FILE * pReadCount = fopen(readCountFile.c_str(), "w");
    fprintf(pReadCount, "%lu\n", readCount/2);
    fclose(pReadCount);
}
#else
    FILE * pReadCount = fopen(readCountFile.c_str(), "w");
    fprintf(pReadCount, "%lu\n", readCount/2);
    fclose(pReadCount);
#endif

#ifdef MPI_ENABLED
    MPI_Finalize();
#endif

    return 0;

}

