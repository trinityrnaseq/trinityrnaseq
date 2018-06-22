#include "sequenceUtil.hpp"
#include <stdlib.h>
#include "stacktrace.hpp"
#include <math.h>
#include <iostream>
#include <sstream>

static const int MAX_LINE_LENGTH = 100000;

char _int_to_base [4] = {'G', 'A', 'T', 'C'};
unsigned char _base_to_int [256] = {
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  20
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  40
		255, 255, 255, 255, 255,   1, 255,   3, 255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, //  60
		255, 255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   1, 255,   3, //  80
		255, 255, 255,   0, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   2, 255, 255, 255, // 100
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220
		255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240
};
//        0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19

string currAccession;
 
bool contains_non_gatc (string kmer) {

  for (unsigned int i = 0; i < kmer.size(); i++) {
	unsigned char c = kmer[i];

	if (_base_to_int[c] > 3)
	  return(true);
/*
	if (! (c == 'g' || c == 'G'
                   || c == 'a' || c == 'A'
                   || c == 't' || c == 'T'
		   || c == 'c' || c == 'C')
		) {
	  return(true);
	}
*/
  }

  return(false);
}


string read_sequence_from_file (string filename) {
  
  ifstream fileReader (filename.c_str());
  if (!fileReader) { // couldn't open file
	throw(stacktrace() + "\n\nCould not open " + filename + "\n");
  }
  
  string mySequence;
  
  char c_line [MAX_LINE_LENGTH];
  
  
  fileReader.getline(c_line, MAX_LINE_LENGTH);
  while (! fileReader.eof() ) {
	if (c_line[0] != '>') {
	  string s_line (c_line);
	  
	  for (unsigned int i = 0; i < strlen(c_line); i++) {
		char c = tolower(c_line[i]);
		if (c == ' ' || c == '\t' || c == '\n') { continue;}
		mySequence += c;
	  }
	}
	
	fileReader.getline(c_line, MAX_LINE_LENGTH);
  }
  
  return (mySequence);
}


fastaRecord readNextFastaRecord (ifstream& reader) {
  
  fastaRecord fr;

  if (reader.eof()) {
	return (fr);
  }
  
  string line;
  
  // prompt till the accession
  if (currAccession.length() == 0) {
	// must read next accession
	while (getline(reader, line)) {
	  if (line[0] == '>') {
		currAccession = line;
		break;
	  }
	}
  }
  
  fr.accession = currAccession;
  
  while (getline(reader, line)) {
	if (line[0] == '>') {
	  currAccession = line;
	  break;
	}
	else {
	  // append characters that are not whitespace
	  for (unsigned int i = 0; i < line.size(); i++) {
		char c = line[i];
                
		if (c != ' ' && c != '\t' && c != '\n') {
		  fr.sequence += c;
		}
	  }
	}
  }
   
  return(fr);
}
        
string revcomp (const string kmer) {
  
  string revstring;
  
  for (int i = kmer.size() -1; i >= 0; i--) {
	char c = kmer[i];
	char revchar;
	
	switch (c) {
	
	case 'g':
	  revchar = 'c';
	  break;
	
	case 'G':
	  revchar = 'C';
	  break;
	  
	case 'a':
	  revchar = 't';
	  break;

	case 'A':
	  revchar = 'T';
	  break;

	case 't':
	  revchar = 'a';
	  break;
	  
	case 'T':
	  revchar = 'A';
	  break;

	case 'c':
	  revchar = 'g';
	  break;
	  
	case 'C':
	  revchar = 'G';
	  break;

	default:
	  revchar = 'N';
	}
	
        
	revstring += revchar;
        
  }
  
  
  return (revstring);
}

kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length)
{
	kmer_int_type_t rev_kmer = 0;
	kmer = ~kmer;
	for (unsigned int i = 0; i < kmer_length; i++) {

		int base = kmer & 3;
		rev_kmer = rev_kmer << 2;
		rev_kmer+= base;
		kmer = kmer >> 2;

	}

	return rev_kmer;
}

kmer_int_type_t get_canonical_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length) {
    // same as DS kmer val
    return(get_DS_kmer_val(kmer_val, kmer_length));
}




string remove_whitespace (string s) {

  string r = "";

  for (unsigned int i = 0; i < s.length(); i++) {
	
	char c = s[i];
	
	if (c != '\t' && c != '\n' && c != ' ') {
	  r += c;
	}
  }

  return(r);
}

int base_to_int_value (char nucleotide) {

  switch (nucleotide) {
	
  case 'G':
  case 'g':
	return(0);
  
  case 'A':
  case 'a':
	return(1);

  case 'T':
  case 't':
	return(2);
	
  case 'C':
  case 'c':
	return(3);
	
  default:
	return(-1);
  
  }

}

char int_to_base(int baseval) {
  
  if (baseval < 0 || baseval > 3) {
	throw (stacktrace() + "\n\nError, baseval out of range 0-3");
  }

  return(_int_to_base[baseval]);
}


kmer_int_type_t kmer_to_intval(string kmer) {
  
  if (kmer.length() > 32) {
	throw(stacktrace() + "\n\nerror, kmer length exceeds 32");
  }

  kmer_int_type_t kmer_val = 0;
  
  for (unsigned int i = 0; i < kmer.length(); i++) {
	unsigned char c = kmer[i];
	int val = _base_to_int[c];
        /* ottmi: Don't need this here: the only possible source for non-gatc characters is
           the reads file and we already check in add_kmer() */
        /* bhaas-to-ottmi: I need to put this back in....  There are other sources where we
           don't explicitly check, and not doing so here causes problems downstream */
        /* deccles-to-bhass: In the interest of code speed / efficiency, the non-gatc check
           is better inlined here -- it's already working out the value per-char, so we get
           the check for a single if statement per character */
        if(val > 3){
          stringstream errstr;
          errstr << "\n\nerror, kmer contains nongatc: " << kmer;
          cerr << errstr.str();
          throw(stacktrace() + errstr.str());
        }
   
	// cerr << "Char " << c << "=" << val << endl;

	kmer_val = kmer_val << 2;
	kmer_val |= val;
	
	//cerr << "\tkmerval: " << kmer_val << endl;
  
	
  }
  
  // cerr << "Kmer: " << kmer << " = " << kmer_val << endl;
  
  return(kmer_val);
}

string decode_kmer_from_intval(kmer_int_type_t intval, unsigned int kmer_length) {

  string kmer(kmer_length, ' ');
  
  for (unsigned int i = 1; i <= kmer_length; i++) {
	
	int base_num = intval & 3ll;
	
	kmer[kmer_length-i] = _int_to_base[base_num];
	
	// cerr << "base: " << base << endl;
	
	intval = intval >> 2;
  }

  return(kmer);
}

float compute_entropy(kmer_int_type_t kmer, unsigned int kmer_length) {

  char counts[] = { 0, 0, 0, 0 };

  for (unsigned int i = 0; i < kmer_length; i++) {

	int c = kmer & 3;
	kmer = kmer >> 2;
	counts[c]++;
  }

  float entropy = 0;

  for (unsigned int i = 0; i < 4; i++) {

	float prob = (float)counts[i] / kmer_length;

	if (prob > 0) {
	  float val = prob * log(1/prob)/log(2.0f);
	  entropy += val;
	}
  }

  return(entropy);
}
  

float compute_entropy(string& kmer) {

  map<char,int> char_map;

  for (unsigned int i = 0; i < kmer.length(); i++) {
	
	char c = kmer[i];
	char_map[c]++;
  }

  float entropy = 0;
  
  char nucs[] = { 'G', 'A', 'T', 'C' };

  for (unsigned int i = 0; i < 4; i++) {
	
	char nuc = nucs[i];
	
	int count = char_map[nuc];
	
	float prob = (float)count / kmer.length();
	
	if (prob > 0) {
	  float val = prob * log(1/prob)/log(2.0f);
	  entropy += val;
	}
  }

  return(entropy);
}



kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length) {
    
    kmer_int_type_t rev_kmer = revcomp_val(kmer_val, kmer_length);
    
    if (rev_kmer > kmer_val)
        kmer_val = rev_kmer;
    
    return(kmer_val);
    
}

vector<kmer_int_type_t> sequence_string_to_kmer_int_type_vector(const string& sequence, int kmer_length) {

    vector<kmer_int_type_t> kit_vec;
    
    for (unsigned int i = 0; i <= sequence.length() - kmer_length; i++) {
        
        string kmer = sequence.substr(i, kmer_length);
                
        kmer_int_type_t k = kmer_to_intval(kmer);
        
        kit_vec.push_back(k);
    }

    return(kit_vec);
    
}

string replace_nonGATC_chars_with_A (string& str) {
    
    stringstream newStr;

    for (unsigned int i = 0; i < str.size(); i++) {
        unsigned char c = str[i];
        
        if (_base_to_int[c] > 3)
            c = 'A';
        
        newStr << c;
    }

    return(newStr.str());

}


unsigned long long generateHash(string& s) {

    // adapted from: http://stackoverflow.com/questions/8094790/how-to-get-hash-code-of-a-string-in-c

    unsigned long long combined_hashcode = 0;

    unsigned int hash = 0;
    for(size_t i = 0; i < s.length(); ++i) {
        char nucleotide = s[i];
        hash = 65599 * hash + nucleotide;

        int base_val = base_to_int_value(nucleotide) + 1; // adding one just in case, since non-gatc = -1, but shouldn't encounter non-gatc here anyway.
        combined_hashcode += base_val;
        
    }
    hash = hash ^ (hash >> 16);

    combined_hashcode <<= 32;
    
    combined_hashcode |= hash;

    return(combined_hashcode);
    
}
