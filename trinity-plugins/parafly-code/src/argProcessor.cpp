#include "argProcessor.hpp"
#include <stdlib.h>
// #include <cstdlib>

ArgProcessor::ArgProcessor(int argc, char* argv[]) {
  
  for (int i=1; i < argc; i++) {
	char* arg = argv[i];
	string myArg (arg);
	if (arg[0] == '-') { 
	
	  argSet[myArg] = true;
  	  	  
	  // see if next argument is another arg or value:
	  if (i != argc-1) {
		// value:
		string nextArg (argv[i+1]);
		argVal[myArg] = nextArg;
	  }
	}
  }
  
}


bool ArgProcessor::isArgSet(string arg) {
  
  map<string,bool>::const_iterator it;

  it = argSet.find(arg);
  if (it != argSet.end()) {
	return(true);
  }
  else {
	return(false);
  }
}

int ArgProcessor::getIntVal(string arg) {
  
  return(atoi(argVal[arg].c_str()));
  
}

long ArgProcessor::getLongVal(string arg) {
  return(atol(argVal[arg].c_str()));
}


float ArgProcessor::getFloatVal(string arg) {
  
  return(atof(argVal[arg].c_str()));
}


string ArgProcessor::getStringVal(string arg) {
  
  return(argVal[arg]);
}


