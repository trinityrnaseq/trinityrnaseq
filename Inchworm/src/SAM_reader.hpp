#ifndef __SAM_READER__

#define __SAM_READER__

#include <fstream>
#include <iostream>
#include <string>
#include "SAM_entry.hpp"

using namespace std;

class SAM_reader {


public:
  SAM_reader(string sam_filename);

  SAM_entry get_next();

  bool has_next();

private:
  
  ifstream filereader;  
  string next_line;

  void advance_filereader();

};

#endif


  
