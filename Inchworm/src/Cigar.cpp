#include "Cigar.hpp"
#include <string>
#include <sstream>
#include "sequenceUtil.hpp"
#include <stdlib.h>

namespace Cigar {

  vector<t_cigar> parse_cigar(string& cigar_txt) {
	
	vector<t_cigar> cigar_eles;
	
	string len_text;
	
	for (unsigned int i = 0; i < cigar_txt.length(); i++) {
	  char c = cigar_txt[i];
	  if (c >= '0' && c <= '9') {
		len_text = len_text + c;
	  }
	  else {
		t_cigar ele;
		ele.len = atol(len_text.c_str());
		ele.code = c;
		cigar_eles.push_back(ele);
		len_text = "";
	  }
	}
    
	return(cigar_eles);
	
  }



  string construct_cigar(vector<alignment_segment>& segments, unsigned int read_length, 
						 const string& genome_seq, char strand) {

	stringstream cigar_text;

	bool check_splice = (genome_seq != "" && strand);
	
	for (unsigned int i = 0; i < segments.size(); i++) {
	  
	  alignment_segment seg = segments[i];
	  
	  unsigned long curr_genome_lend = seg.genome_lend;
	  unsigned long curr_genome_rend = seg.genome_rend;
	  unsigned long curr_query_lend = seg.query_lend;
	  unsigned long curr_query_rend = seg.query_rend;
	  
	  if (i == 0) {

		if (curr_query_lend > 1) {
		  
		  cigar_text << curr_query_lend - 1 << "S";
		}
	  }
	  else {

		alignment_segment prev_seg = segments[i-1];
		
		unsigned long prev_genome_rend = prev_seg.genome_rend;
		unsigned long prev_query_rend = prev_seg.query_rend;

		unsigned long delta_genome = curr_genome_lend - prev_genome_rend;
		if (delta_genome > 1) {
		  char deletion_intron_char = (check_splice) 
			? check_intron_consensus(prev_genome_rend, curr_genome_lend, genome_seq, strand)
			: 'D';

		  cigar_text << delta_genome - 1 << deletion_intron_char;
		}

		unsigned long delta_query = curr_query_lend - prev_query_rend;
		if (delta_query > 1) {
		  cigar_text << delta_query - 1 << 'I';
		}
	  }
		
	  unsigned long len = curr_genome_rend - curr_genome_lend + 1;
		
	  cigar_text << len << 'M';
	  
	  if (i == segments.size() - 1) {
		
		if (curr_query_rend < read_length) {
		  cigar_text << read_length - curr_query_rend << 'S';
		}
	  }
	  
	}
	
	return(cigar_text.str());
	
  }


  char check_intron_consensus(unsigned long prev_exon_rend, unsigned long next_exon_lend, 
							  const string& genome_seq, char strand) {

	
	unsigned int genome_seq_length = genome_seq.length();
	if (prev_exon_rend > genome_seq_length || next_exon_lend > genome_seq_length) {
	  stringstream errmsg;
	  errmsg << "Error, coordinates: " << prev_exon_rend << " and " << next_exon_lend 
			 << " are not entirely within genome sequence length: " << genome_seq_length;

	  throw(errmsg.str());
	}
	
	
	string left_dinuc;
	left_dinuc += genome_seq[prev_exon_rend];
	left_dinuc += genome_seq[prev_exon_rend + 1];
	
	string right_dinuc;
	right_dinuc += genome_seq[next_exon_lend -3];
	right_dinuc += genome_seq[next_exon_lend-2];
	
	if (strand == '-') {
	  string left_dinuc_copy = left_dinuc;
	  string right_dinuc_copy = right_dinuc;
	  
	  left_dinuc = revcomp(right_dinuc_copy);
	  right_dinuc = revcomp(left_dinuc_copy);
	}

	if (
		( (left_dinuc == "GT" || left_dinuc == "GC") && right_dinuc == "AG")
		||
		(left_dinuc == "CT" && right_dinuc == "AC")
		) {
	  
	  return('N'); // has proper splice boundaries.
	}
	else {
	  return('D');
	}

  }
  
  


}  // end of namespace 


