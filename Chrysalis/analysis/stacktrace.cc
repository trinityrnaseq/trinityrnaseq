#include "analysis/stacktrace.h"
#include <signal.h>

#include <sstream>
// #include <cxxabi.h>
#include <vector>
#include <stdlib.h>


#ifdef HAVE_STACKTRACE_H

#include <execinfo.h>


// private prototypes
void _tokenize(const string& str,
			   vector<string>& tokens,
			   const string& delimiters = " ");

// implementation

string stacktrace() {
  
  /* based on article described here:
	 http://www.linuxjournal.com/article/6391
	 Stack Backtracing Inside Your Program
	 Aug 11, 2003  By Gianluca Insolvibile
	 Linux Journal
  */
  
  void *trace[16];
  char **messages = (char **)NULL;
  int i, trace_size = 0;
  
  // int status;
  
  stringstream s;
  
  trace_size = backtrace(trace, 16);
  messages = backtrace_symbols(trace, trace_size);
  
  s << endl << "[bt] Execution path:" << endl;
  
  for (i=0; i<trace_size; ++i) {
    string message = messages[i];
		
	// construct backtrace output

    s << "[bt] " << message << endl;

  }
  
  return(s.str());
}


#else

string stacktrace () {
  return(" <stacktrace off> ");
}

#endif



