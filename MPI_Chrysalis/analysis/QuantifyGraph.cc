//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include <string.h>

#include "base/CommandLineParser.h"

#include "aligns/KmerAlignCore.h"
#include "analysis/DNAVector.h"
#include <math.h>
#include "base/FileParser.h"
#include "analysis/KmerTable.h"

#include "analysis/CompMgr.h"

extern "C"
{
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
}


static bool DEBUG = false;


// need to add Execute() and Exists() to a common library, since used in multiple chrysalis components (including Chrysalis.cc)
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



void SortPrint(FILE * pReads, svec<IDS> & ids, const vecDNAVector & seq) 
{
    long long i;
    
    Sort(ids);
    int lastID = -1;
    int id = -1;
    int start = -1;
    int edge = -1;
    int lastStart = -1;
    int lastEdge = -1;
    int ori;
    
    string line;
    char tmp[1024];
    
    int lastStartTemp = -1;
    int lastOri = 1;
    
    for (i=0; i<ids.lsize(); i++) {
        id = ids[i].ID();
        ori = ids[i].Ori();
        start = ids[i].Start();
        edge = ids[i].Edge();
        //cout << id << "\t" << start << "\t" << edge << endl;
        if (id != lastID
#ifndef NO_REVERSE_OUT	  
            || ori != lastOri
#endif
            ) {
            if (lastID != -1) {
                
                sprintf(tmp, "%d\t%d\t", lastStart, lastEdge);	
                line += tmp;
                if (lastStart > lastStartTemp) {
                    fprintf(pReads, "%s\t", line.c_str());
                    //const DNAVector &d = seq[lastID];
                    
#ifndef NO_REVERSE_OUT	  
                    DNAVector d = seq[lastID];
                    if (lastOri == -1) {
                        d.ReverseComplement();
                        //cout << "Reversing" << endl;
                    } else {
                        //cout << "Forward" << endl;
                    }
#endif 
                    
                    for (int j=0; j<d.isize(); j++) {
                        tmp[1] = 0;
                        tmp[0] = d[j];
                        fprintf(pReads, "%s", tmp);
                    }
                    
                    if (lastOri == -1)
                        fprintf(pReads, "\t-");
                    else
                        fprintf(pReads, "\t+");
                    fprintf(pReads, "\n");
                }
                //fprintf(pReads, "%d\t%d\n", lastStart, lastEdge);
                line = "";
            }
            //fprintf(pReads, "%s\t%d\t%d\t", seq.Name(id).c_str(), start, edge);
            sprintf(tmp, "%s\t%d\t%d\t", seq.Name(id).c_str(), start, edge);
            line = tmp;
            
            lastStartTemp = start;
        }
        lastID = id;
        lastStart = start;
        lastEdge = edge;
        lastOri = ori;
        
    }
    
    if (id != -1) {
        sprintf(tmp, "%d\t%d\t", start, edge);	
        line += tmp;
        if (lastStart > lastStartTemp) {
            fprintf(pReads, "%s\t", line.c_str());
            DNAVector d = seq[id];
            if (ori == -1)
                d.ReverseComplement();
            for (int j=0; j<d.isize(); j++) {
                tmp[1] = 0;
                tmp[0] = d[j];
                fprintf(pReads, "%s", tmp);
            }
            if (lastOri == 1)
                fprintf(pReads, "\t+");
            else
                fprintf(pReads, "\t-");
            
            
            fprintf(pReads, "\n");
        }
        //fprintf(pReads, "%d\t%d\n", start, edge);
    }
}

//========================================================================
//========================================================================
//========================================================================


bool Irregular(char l)
{
    if (l == 'A' || l == 'C' || l == 'G' || l == 'T')
        return false;
    //cout << "Irregular char: " << l << endl;
    return true;
}



string ReadsExt(const string & in) 
{
    char tmp[1024];
    strcpy(tmp, in.c_str());
    int n = strlen(tmp);
    
    
    for (int i=n-1; i>=0; i--) {
        if (n-i > 6) {
            break;
        }
        if (tmp[i] == '.') {
            tmp[i] = 0;
            string out = tmp;
            out += ".reads";
            return out;
        }
        
    }
    string out = in + ".reads";
    return out;
}

int main(int argc,char** argv)
{
    
    commandArg<string> aStringCmmd("-i","read fasta file");
    commandArg<string> gStringCmmd("-g","graph file");
    commandArg<string> oStringCmmd("-o","graph output");
    commandArg<int> kCmmd("-k","kmer size", 24);
    commandArg<bool> strandCmmd("-strand","strand specific", false);
    //commandArg<bool> cCmmd("-nc","do not fully connected graph", false);
    commandArg<long> maxReadsCmd("-max_reads", "max number of reads to map to graph", -1);
    commandArg<bool> debugCmmd("-debug", "verbosely describe operations", false);
    commandArg<bool> no_cleanupCmmd ("-no_cleanup", "retain input files on success", false);
    
    commandLineParser P(argc,argv);
    P.SetDescription("Assembles k-mer sequences.");
    P.registerArg(aStringCmmd);
    P.registerArg(gStringCmmd);
    P.registerArg(oStringCmmd);
    P.registerArg(kCmmd);
    P.registerArg(strandCmmd);
    //P.registerArg(cCmmd);
    P.registerArg(maxReadsCmd);
    P.registerArg(debugCmmd);
    P.registerArg(no_cleanupCmmd);

    P.parse();
    
    string aString = P.GetStringValueFor(aStringCmmd); // reads
    string gString = P.GetStringValueFor(gStringCmmd); // graph input
    string oString = P.GetStringValueFor(oStringCmmd); // graph output
    bool sStrand = P.GetBoolValueFor(strandCmmd);
    int k = P.GetIntValueFor(kCmmd)+1;
    long max_reads = P.GetLongValueFor(maxReadsCmd); 
    bool NO_CLEANUP = P.GetBoolValueFor(no_cleanupCmmd);
    
    DEBUG = P.GetBoolValueFor(debugCmmd);
    


    
    if (Exists(oString) && (! Exists(gString)) && (! Exists(aString)) ) {
        cerr << "Quantify graph previously finished successfully on " << aString << ".  Not rerunning here." << endl;
        return(0);
    }
    else if (! (Exists(gString) && Exists(aString)) ) {
        cerr << "ERROR: missing either: " << gString << " or " << aString << ", cannot run QuantifyGraph here." << endl;
        return(1);
    }
    
    
    int i, j;
    //vecbasevector contigBases;
    
    vecDNAVector seq;

    if (max_reads > 0) {
        // std::cerr << "*Restricting number of input reads to " << max_reads << endl;
        seq.setMaxSeqsToRead(max_reads);
    }
    
    seq.Read(aString, false, true, true, 1000); // parse the reads from the fasta file
    
    KmerSequence kmers(k, &seq); // kmer catalog based on the fasta reads
    kmers.Add(seq);
    
    long long m = kmers.GetBoundValue();
    
    FlatFileParser parser; // read the raw graph
    parser.Open(gString);
    
    FILE * pOut = fopen(oString.c_str(), "w");  // output graph
    string reads = ReadsExt(oString);
    FILE * pReads = fopen(reads.c_str(), "w");  // output reads in context of graph
    
    svec<IDS> ids;
    ids.reserve(seq.isize());
    
    //string last;
    //int lastNode = -1;
    
    svec<char> first;
    first.resize(100000, 'N');

    // do an initial scan to set up the node identities and linkage info
    while (parser.ParseLine()) {
        
        if (parser.GetItemCount() >= 4) {
            const string & s = parser.AsString(3); // kmer
            int node = parser.AsInt(0);
            int prevNode = parser.AsInt(1);
            
            const char * p2 = s.c_str();
            if (node >= first.isize())
                first.resize(node + 10000, 'N');

            first[node] = p2[0]; // first letter of the kmer stored
        }
    }

    
    // now, do a second pass:
    parser.Open(gString);

    while (parser.ParseLine()) {
        
        if (parser.GetItemCount() < 4) {

            fprintf(pOut, "%s\n", parser.Line().c_str()); // component header line

            // processing of component data from previously processed component
            
            if (ids.lsize() > 0) {
                SortPrint(pReads, ids, seq);
                //UniqueSort(ids);
                //for (i=0; i<ids.isize(); i++) {
                //fprintf(pReads, "%s\n", seq.Name(ids[i]).c_str());
                //}
            }
            fprintf(pReads, "%s\n", parser.Line().c_str());
            ids.clear();
            
            //first.clear();
            //first.reserve(50000);
            continue;
        }
        
        const string & s = parser.AsString(3); // kmer
        int node = parser.AsInt(0);
        int prevNode = parser.AsInt(1);
        
        const char * p2 = s.c_str();
        //if (node >= first.isize())
        //    first.resize(node + 10000, 'N');
        //first[node] = p2[0];
        
        //if (prevNode >= 0 && last == "" ) {
        //  cout << "Potential ERROR" << endl;
        //}
        //if(prevNode >= 0 && prevNode != lastNode) {
        //  cout << "Potential ERROR (2): prevNode = " << prevNode;
        //  cout << " lastNode=" << lastNode << endl;      
        //}
        //int edge = parser.AsInt(0);
        long long edge = prevNode;
        long long n1 = 0;
        long long n2 = 0;
        if (prevNode >= 0) {
            
            //  building the whole kmer sequence in 'sub'
            DNAVector sub;
            sub.resize(strlen(s.c_str())+1);
            const char * p = s.c_str();
            for (i=0; i<sub.isize()-1; i++)
                sub[i+1] = p[i];
            
            //const char * p2 = last.c_str();
            if (first[prevNode] == 'N')
                cout << "ERROR!! first[prevNode] where prevNode = " << prevNode << " unset" << endl;

            sub[0] = first[prevNode];
            


            kmers.BasesToNumberCountPlus(ids, n1, sub, edge);
            
            if (!sStrand) {
                sub.ReverseComplement();
                long long from = ids.lsize();
                kmers.BasesToNumberCountPlus(ids, n2, sub, edge);
                
                if (n1 + n2 < 0x7FFFFFFF) {
                    
                    for (long long x=from; x<ids.lsize(); x++) {
                        ids[x].SetOri(-1);
                        int len = seq[ids[x].ID()].isize();
                        int pos = ids[x].Start()+1;
                        //cout << "len=" << len << " pos=" << pos;
#ifndef NO_REVERSE_OUT	  
                        ids[x].SetStart(len-pos-k+1);
                        //cout << " new=" << len-pos-k+1 << endl;
#else
                        ids[x].SetStart(pos+1);
#endif
                    }
                } else {
                    cout << "WARNING: k-mer overflow, n=" << n1 + n2 << ". Discarding." << endl;
                    n1 = n2 = 0;
                    //ids.resize(0);	
                }
            }
            
            /*  if (s == "ATATCACAAAACAATCTTCATTCG") {
                for (int x = 0; x<sub.isize(); x++)
                cout << sub[x];
                cout << endl;
                
                cout << "Count=" << n1 + n2 << endl;
                }*/
        }
        
        for (i=0; i<parser.GetItemCount(); i++) {
            if (i>0)
                fprintf(pOut, "\t");
            if (i == 2) {
                fprintf(pOut, "%d", (int)(n1+n2));
            } else {
                fprintf(pOut, "%s", parser.AsString(i).c_str());
            }
        }
        fprintf(pOut, "\n");
        //lastNode = node;
        //last = s;
    }
    if (ids.lsize() > 0) {
        SortPrint(pReads, ids, seq);
    }
    
    fclose(pOut);
    fclose(pReads);

           
    
    
    // only remove the input files once the outputs have been successfully generated.
    if (! NO_CLEANUP) {

        // remove inputs to reduce file counts.
        unlink(aString.c_str());
        unlink(gString.c_str());
    }

    return 0;
    
}
