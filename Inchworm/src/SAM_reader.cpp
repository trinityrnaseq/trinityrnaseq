#include "SAM_reader.hpp"
#include "SAM_entry.hpp"
#include <sstream>

SAM_reader::SAM_reader(string sam_filename) {
  
  filereader.open(sam_filename.c_str());
  if (! filereader.is_open()) {
	stringstream errstr;
	errstr << "Error, cannot open filename: " << sam_filename;
	throw(errstr.str());
	
  }
  
  advance_filereader();
  
  
}



void SAM_reader::advance_filereader() {

  next_line = "";
  
  while (! filereader.eof()) {
	getline(filereader, next_line);
	if (next_line[0] != '@') {
	  break;
	}
  }

}



SAM_entry SAM_reader::get_next() {

  string ret_line = next_line;
  
  advance_filereader();
  
  SAM_entry s(ret_line);
  
  return(s);
}

bool SAM_reader::has_next() {
  if (next_line == "" && filereader.eof()) {
	return(false);
  }
  else {
	return(true);
  }
}

