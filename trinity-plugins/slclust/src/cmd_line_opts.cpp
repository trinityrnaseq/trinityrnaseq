#include "cmd_line_opts.h"
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "common.h"
#include <cstdio>

using namespace std;

/* 
   The cmd-line processing code was borrowed from RepeatScout with 
   permission from Neil Jones
*/


bool co_get_int(int argc, char** argv, const char* text, int* res) {
  int opt_len = strlen(text);
  bool found = false;
  
  int ii;
  for(ii = 1; ii < argc; ++ii) {
    if (VERBOSE >= debug) {
      fprintf (stderr, "co_get_int(%s), comparing to %s\n", text, argv[ii]);
    }
    
    if( argv[ii][0] == '-' ) {
      if(strncmp(text, argv[ii], opt_len) == 0) {
        if ( ii+1 < argc) {
          *res = atoi(argv[ii+1]);
          found = true;
          break;
        }
      }
      
    }
  }

  if (VERBOSE >= debug) {
    cerr << "found: " << found << endl;
  }
  
  return found;
}

bool co_get_bool(int argc, char** argv, const char* text) {
  int opt_len = strlen(text);
  bool found = false;
  
  int ii;
  for(ii = 1; ii < argc; ++ii) {
    
    if (VERBOSE >= debug) {
      fprintf (stderr, "co_get_bool(%s), comparing to %s\n", text, argv[ii]);
    }
    
    
    if( argv[ii][0] == '-' ) {
      if(strncmp(text, argv[ii], opt_len) == 0) {
        found = true;
        break;
      }
      
    }
  }
  
  if (VERBOSE >= debug) {
    cerr << "found: " << found << endl;
  }
  
  return found;
}

bool co_get_float(int argc, char** argv, const char* text, float* res) {
  int opt_len = strlen(text);
  bool found = false;
  
  int ii;
  for(ii = 1; ii < argc; ++ii) {
    
    if (VERBOSE >= debug) {
      fprintf (stderr, "co_get_float(%s), comparing to %s\n", text, argv[ii]);
    }


    if( argv[ii][0] == '-' ) {
      if(strncmp(text, argv[ii], opt_len) == 0) {
        if (ii+1 < argc) {
          *res = atof(argv[ii+1]);
          found = true;
          break;
        }
      }
      
    }
  }
  
  if (VERBOSE >= debug) {
    cerr << "found: " << found << endl;
  }
  return found;
}


bool co_get_string(int argc, char** argv, const char* text, char** res) {
  int opt_len = strlen(text);
  bool found = false;
  
  int ii;
  for(ii = 1; ii < argc; ++ii) {
    
    if (VERBOSE >= debug) {
      fprintf (stderr, "co_get_string(%s), comparing to %s\n", text, argv[ii]);
    }

    if( argv[ii][0] == '-' ) {
      if(strncmp(text, argv[ii], opt_len) == 0) {
        if (ii+1 < argc) {
          *res = argv[ii+1];
          found = true;
          break;
        }
      }
      
    }
  }
 
  if (VERBOSE >= debug) {
    cerr << "found: " << found << endl;
  } 
  return found;
}
