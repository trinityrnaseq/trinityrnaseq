#include "graph.h"
#include <cstdio>
#include <iostream>
#include "cmd_line_opts.h"
#include "common.h"
#include <cstdlib>

using namespace std;

// verbosity settings.
VERBOSITY VERBOSE = none;

// prototypes
void die_usage(string progname);




int main (int argc, char** argv) {
  
  string progname = argv[0];
  
  // cmd line processing
  
  if (
      co_get_bool(argc, argv, "-h") 
      ||
      co_get_bool(argc, argv, "--h") 
      ) {
    
    die_usage(progname);
    
  }
 
  // check verbosity setting
  if (co_get_bool(argc, argv, "-vv")) {
    VERBOSE = debug;
  }
  else if (co_get_bool(argc, argv, "-v")) {
    VERBOSE = info;
  }
  
  // get the Jaccard coefficient if specified.
  float cutoff;
  bool want_jaccard = co_get_float(argc, argv, "-j", &cutoff);
  
  // Begin Graph Building.
  
  Graph g;
  string a,b;
  char line [10000];
  char achars [200];
  char bchars [200];
  int numread = 0;
  char* readLine;
  
  int count = 0;
  
  if (VERBOSE >= info) 
    cerr << "-reading pairs." << endl;
  
  while ( (readLine = fgets (line, 10000, stdin)) != NULL) {
    
    numread = sscanf (readLine, "%s %s", achars, bchars);
    
    if (VERBOSE >= debug)
      printf ("numread: %d\n", numread);
    
    if (numread == 2) {
      
      if (VERBOSE >= debug)
        printf ("Read pairs: (%s, %s)\n", achars, bchars);
      
      a = achars;
      b = bchars;
      
      if (a.compare(b) != 0) {
        
        count++;
        if (VERBOSE >= debug) { 
          cerr << "\r  pairs read: " << count << "     "; 
        }
        else if (VERBOSE >= info && (count % 1000 == 0) ) {
          cerr << "\r pairs read: " << count << "     ";
        }
        
        g.addLinkedNodes(a,b);
                
      }
      
    } else {
      string sline = line;
      cerr << "\nError, line: (" << sline << ") has unrecognized formatting." << endl;
    }
  }
  
  if (VERBOSE >= info) 
    cerr << "\nDone reading pairs.  Now building clusters." << endl;
  

  
  if (want_jaccard) {
    
    /*   Getting Jaccard Clusters by removing edges below the specified coefficent */
    
    if (VERBOSE >= info) 
      cerr << endl << endl << "Finding jaccard clusters at coeff(" << cutoff << ")" << endl;
    
    Graph* g2 = g.applyJaccardCoeff (cutoff);
    
    if (VERBOSE >= info)
      cerr << "Extracting clusters now." << endl;

    g2->printClusters();
    delete g2;
    
    if (VERBOSE >= info) 
      cerr << "Finished." << endl;
    
  }
  else {
    
    /*   Getting Simple Single-Linkage Clusters     */
    
    g.printClusters();
  }
  
  
  return (0);
  
}


void die_usage (string progname) {
  
  cerr << "usage: " << progname << " [opts] < file_of_pairs > clusters " << endl << endl
       << "options:" << endl
       << "\t -j jaccard_coefficient " << endl
       << "\t -v[v] verbosity at 'info', 'debug'  " << endl
       << endl << endl ;
  
  exit(1);
}

