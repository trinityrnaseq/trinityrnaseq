#ifndef __CIGAR__

#define __CIGAR__

#include <vector>
#include <string>
#include "SAM_entry.hpp"

using namespace std;

namespace Cigar {

  typedef struct { 
	unsigned long len;
	char code;
  } t_cigar;
  
  vector<t_cigar> parse_cigar(string& cigar_txt);
  
  string construct_cigar(vector<alignment_segment>& segments, 
						 unsigned int read_length,
						 const string& genome_seq,
						 char strand = NULL);
  
  
  
  char check_intron_consensus(unsigned long prev_exon_rend, unsigned long next_exon_lend, 
							  const string& genome_seq, char strand);
  
  
  

}


#endif
