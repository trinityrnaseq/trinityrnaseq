#include "Fasta_entry.hpp"
#include "string_util.hpp"


// constructor
Fasta_entry::Fasta_entry(string header, string sequence) {
  
  // piece apart the header
  
  if (header[0] == '>') {
	header.erase(0,1);
  }

  
  vector<string> toks;
  string acc;
  
  string_util::tokenize(header, toks, " \t");
  
  if (toks.size() > 0) {
	acc = toks[0];
  }
  
  this->_header = header;
  this->_accession = acc;
  this->_sequence = sequence;
 

}

string Fasta_entry::get_accession() {
  return(this->_accession);
}

string Fasta_entry::get_header() {
  return(this->_header);
}

string Fasta_entry::get_sequence() {
  return(this->_sequence);
}


  
 
  
