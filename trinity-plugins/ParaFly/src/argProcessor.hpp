#ifndef __argProcessor__
#define __argProcessor__

#include <map>
#include <string>
using namespace std;

class ArgProcessor {
  
public:
  ArgProcessor(int argc, char* argv[]);
  
  bool isArgSet(string arg);
  
  int getIntVal(string arg);
  
  long getLongVal(string arg);

  float getFloatVal(string arg);
  
  string getStringVal(string arg);
  
private:
  map<string,bool> argSet;
  map<string,string> argVal;

};


#endif

