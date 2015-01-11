#include <stdlib.h>
#include "SAM_entry.hpp"
#include "string_util.hpp"
#include <iostream>
#include <sstream>
#include "Cigar.hpp"

using namespace Cigar;

// constructor
SAM_entry::SAM_entry(string sam_line) {

  line = sam_line;

  string_util::tokenize(line, tokens, "\t");



}

string SAM_entry::get_sam_string() {
  
  // reconstruct the line
  stringstream sam_line;
  
  for (unsigned int i = 0; i < tokens.size(); i++) {
	sam_line << tokens[i];
	
	if (i < tokens.size() - 1) {
	  sam_line << "\t";
	}
  }
  
  return(sam_line.str());  
}


string SAM_entry::toString() {

  return(get_sam_string());
}



string SAM_entry::get_read_name() {
  return(tokens[0]);
}


string SAM_entry::get_full_read_name() {

    if (is_paired()) {
        
        string core_read_name = get_read_name();
        if (is_first_in_pair()) {
            return(core_read_name + "/1");
        }
        else {
            // assert is second in pair
            return(core_read_name + "/2"); 
        }
        

    }
    else {
        return(get_read_name());
    }


}


string SAM_entry::get_scaffold_name() {
  return(tokens[2]);
}

unsigned long SAM_entry::get_scaffold_position() {
  
  return(atol(tokens[3].c_str()));
}

void SAM_entry::set_scaffold_position(unsigned long position) {
  
  stringstream s;
  s << position;

  tokens[3] = s.str();
  
}

string SAM_entry::get_cigar_alignment() {
  return(tokens[5]);
}

void SAM_entry::set_cigar_alignment(string cigar_string) {
  tokens[5] = cigar_string;
}


vector<alignment_segment> SAM_entry::get_alignment_coords() {

  vector<alignment_segment> segs;

  unsigned long genome_lend = get_scaffold_position();

  string alignment = get_cigar_alignment();

  unsigned long query_lend = 0;

  genome_lend--; // move pointer just before first position

  vector<t_cigar> cigar_structs = parse_cigar(alignment);
  
  for (unsigned int i=0; i < cigar_structs.size(); i++) {

	t_cigar cigar_struct = cigar_structs[i];
	
	unsigned long len = cigar_struct.len;
	char code = cigar_struct.code;

	

	if (code == 'M') {
	  
	  alignment_segment s;
	  
	  unsigned long genome_rend = genome_lend + len;
	  unsigned long query_rend = query_lend + len;
	  
	  s.genome_lend = genome_lend + 1;
	  s.genome_rend = genome_rend;
	  s.query_lend = query_lend + 1;
	  s.query_rend = query_rend;

	  segs.push_back(s);

	  genome_lend = genome_rend;
	  query_lend = query_rend;
	  
	}
	else if (code == 'D' || code == 'N') {
	  genome_lend += len;
	  
	}
	else if (code == 'I' || code == 'S' || code == 'H') {
	  query_lend += len;
	}
	else {

	  stringstream errmsg;
	  errmsg << "Error, cannot parse cigar code: " << code;
	  throw(errmsg.str());
	}
	
  }
  
  return(segs);
}

string SAM_entry::get_mate_scaffold_name() {
  return(tokens[6]);
}

unsigned long SAM_entry::get_mate_scaffold_position() {
  return(atol(tokens[7].c_str()));
}

unsigned int SAM_entry::get_mapping_quality() {
  return(atoi(tokens[4].c_str()));
}

string SAM_entry::get_sequence() {
  return(tokens[9]);
}

string SAM_entry::get_quality_scores() {
  return(tokens[10]);
}

//***** FLAG Processing ******//

unsigned int SAM_entry::get_flag() {
  return(atoi(tokens[1].c_str()));
}

bool SAM_entry::is_paired() {
  return(get_flag() & 0x0001);
}

bool SAM_entry::is_proper_pair() {
  return(get_flag() & 0x0002);
}

bool SAM_entry::is_query_unmapped() {
  return(get_flag() & 0x0004);
}

bool SAM_entry::is_mate_unmapped() {
  return(get_flag() & 0x0008);
}

char SAM_entry::get_query_strand() {
  if (get_flag() & 0x0010) {
	return('-');
  }
  else {
	return('+');
  }
}

char SAM_entry::get_query_transcribed_strand() {
  
  char strand = get_query_strand();
  
  if (is_paired() && is_first_in_pair()) {

	char transcribed_strand = (strand == '+') ? '-' : '+';

	return(transcribed_strand);
  }
  else {
	return (strand);
  }
  
  
  
}

char SAM_entry::get_mate_strand() {
  if (get_flag() & 0x0020) {
	return('-');
  }
  else {
	return('+');
  }
}

bool SAM_entry::is_first_in_pair() {
  return(get_flag() & 0x0040);
}

bool SAM_entry::is_second_in_pair() {
  return(get_flag() & 0x0080);
}


  
