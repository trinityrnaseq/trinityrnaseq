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
#include "analysis/DNAVector.h"



typedef vector<string> string_vec;

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


void write_iworm_bundle (string filename, vector<string_vec>& bundled, vector<int>& component_ids) {

    ofstream ofh;
    ofh.open(filename.c_str());
    

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
        
    }

    ofh.close();

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

    commandArg<string> iStringCmmd("-i","input file: ie. GraphFromIwormFasta.out");
    commandArg<string> oStringCmmd("-o","output file: ie. bundled_iworm_contigs.fasta");
    commandArg<int> minCmmd("-min","minimum sequence length", 300);
    commandArg<bool> debugCmmd("-debug", "verbosely describes operations", false);
    
    commandLineParser P(argc,argv);
    P.SetDescription("Assemble transcriptomes from reads.");
    P.registerArg(iStringCmmd);
    P.registerArg(oStringCmmd);
    P.registerArg(minCmmd);
    P.registerArg(debugCmmd);

    P.parse();

    
    string components_file = P.GetStringValueFor(iStringCmmd);
    string bundledName = P.GetStringValueFor(oStringCmmd);
    int minLen = P.GetIntValueFor(minCmmd);

    
    bool DEBUG = P.GetBoolValueFor(debugCmmd);

    
    ///////////////////////////////////////////////////////////////////////////////////////
    // Read components.out, create chrysalis/bundled.fasta which represents clustered inchworm sequences
    
    
    FlatFileParser parser;  
    parser.Open(components_file);  // read components.out file
    
    int component_no = 0;
    
    //ComponentFileMgr mgr(outDir);
    
    vecDNAVector tmpSeq;
    vector<int> component_ids;
    
    vector<string_vec> bundled;
    //bundled.reserve(1000000);
    
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
            
            
            
            vector<string> iworm_bundle;
            int sum_len = 0;
            for (int x = 0; x < tmpSeq.isize(); x++) {
                iworm_bundle.push_back(get_seq_string(tmpSeq[x]));    
                sum_len += tmpSeq[x].isize();
            }
            
            if (sum_len < minLen) {
                // skipping, too short cumulative seq
                continue;
            }
            
            bundled.push_back(iworm_bundle);
            component_ids.push_back(component_no);
            
            
            tmpSeq.resize(0);
            
        }
        
    }
    
    write_iworm_bundle(bundledName, bundled, component_ids);
        
    
    return(0);

}
