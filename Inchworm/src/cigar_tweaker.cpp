#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "string_util.hpp"
#include <vector>
#include "SAM_reader.hpp"
#include "SAM_entry.hpp"
#include "Cigar.hpp"
#include <map>
#include "Fasta_entry.hpp"
#include "Fasta_reader.hpp"


using namespace std;


unsigned long compute_segment_length(alignment_segment& seg);

string describe_segments(vector<alignment_segment>& v);

void populate_genome_map(string& genome_fa, map<string,string>& genome_map);

int run (int argc, char* argv[]);


int main (int argc, char* argv[]) {

  stringstream usage;
  usage << "usage: " << argv[0] << " file.sam genome.fa [min_terminal_seg_length=10]" << endl << endl;
  
  if (argc < 3) {
	cerr << usage.str();
	exit(1);
	}

  
  try {
	return(run(argc, argv));
  }
  
  catch (string errmsg) {
	
	cerr << "Fatal error: " << errmsg << endl;
	exit(2);
	
  }
  

}

int run (int argc, char* argv []) {

  unsigned long MIN_TERMINAL_SEGMENT_LENGTH = 10;
  if (argc >= 4) {
	MIN_TERMINAL_SEGMENT_LENGTH = atol(argv[3]);
  }
  
  string sam_filename = argv[1];
  string genome_fa = argv[2];
  
  map<string,string> genome_map;
  
  SAM_reader sr(sam_filename);
  
  populate_genome_map(genome_fa, genome_map);
  
  while (sr.has_next()) {
	
	SAM_entry se = sr.get_next();

	string line = se.get_sam_string();
	
	string genome_acc = se.get_scaffold_name();
	
	string sequence = se.get_sequence();


	string original_cigar = se.get_cigar_alignment();

	stringstream report;

	if (! se.is_query_unmapped()) {
	  vector<alignment_segment> v = se.get_alignment_coords();
	  
	  // check to see if we need to update intron information or adjust boundaries.
	  
	  

	  if (v.size() > 1) {
		
		report << endl << original_cigar << endl;
		report << describe_segments(v);
		vector<alignment_segment> new_v;
		
		// first do a terminal segment length check
		bool exclude_first = (compute_segment_length(v[0]) < MIN_TERMINAL_SEGMENT_LENGTH);
		bool exclude_last = (compute_segment_length(v[v.size()-1]) < MIN_TERMINAL_SEGMENT_LENGTH);
		
		if (exclude_first || exclude_last) {
		  for (unsigned int i = 0; i < v.size(); i++) {
			if (i == 0 && exclude_first) { 
			  continue;
			}
			if (i == v.size()-1 && exclude_last) {
			  continue;
			}
			
			new_v.push_back(v[i]);
		  }
		
		  // update cigar string.
		  string new_cigar = Cigar::construct_cigar(new_v, sequence.length(), "");
		  se.set_cigar_alignment(new_cigar);
		  if (exclude_first) {
			unsigned long new_scaff_start = new_v[0].genome_lend;
			se.set_scaffold_position(new_scaff_start);
		  }
		}
		else {
		  new_v = v;
		}
		
		// run splice check
		
		if (new_v.size() > 1) {


		  if (genome_map.find(genome_acc) == genome_map.end()) {
			stringstream errmsg;
			errmsg << "Error, no genome sequence available for: " << genome_acc << endl;
			throw(errmsg.str());
		  }

		  const string& genome = genome_map[ genome_acc ];
		  
		  string new_cigar = Cigar::construct_cigar(new_v, sequence.length(), genome, se.get_query_transcribed_strand());
		  
		  se.set_cigar_alignment(new_cigar);
		}
		
	  } // end of if size > 1
	  
	

	  // report SAM line.
	  cout << se.toString() << endl;
	
	
	}
  }	
  return(0);
}

  
unsigned long compute_segment_length(alignment_segment& seg)  {

  return(seg.genome_rend - seg.genome_lend + 1);

}

	
string describe_segments(vector<alignment_segment>& v) {

  stringstream text;

  for (unsigned int i = 0; i < v.size(); i++) {
	
	alignment_segment s = v[i];
	
	text << s.genome_lend << "-" << s.genome_rend << "\t"
		 << s.query_lend << "-" << s.query_rend << endl;
  }

  return(text.str());

}


void populate_genome_map(string& genome_fa, map<string,string>& genome_map) {

  Fasta_reader fr(genome_fa);

  while (fr.hasNext()) {
	Fasta_entry fe = fr.getNext();

	string accession = fe.get_accession();
	// cerr << "storing: [" << accession << "]" << endl;
	
	genome_map[ accession ] = fe.get_sequence();
	
  }

}

