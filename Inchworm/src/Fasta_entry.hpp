#ifndef __FASTA_ENTRY__
#define __FASTA_ENTRY__

#include <string>

using namespace std;

class Fasta_entry {

 public:

  Fasta_entry() {};
  Fasta_entry(string header, string sequence);

  string get_accession();

  string get_header();

  string get_sequence();

  
 private:
  string _accession;
  string _header;
  string _sequence;
  
};

#endif

