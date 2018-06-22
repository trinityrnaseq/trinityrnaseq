#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <map>
#include <cctype>



#ifndef __SEQUENCEUTIL__

#define __SEQUENCEUTIL__

using namespace std;

// misc typedefs
typedef struct { string accession; string header; string sequence; } fastaRecord;

typedef unsigned long long kmer_int_type_t;

// function prototypes
string read_sequence_from_file (string filename);
string revcomp (const string);
fastaRecord readNextFastaRecord(ifstream& reader);
bool contains_non_gatc(string kmer);
string remove_whitespace(string s);


char int_to_base(int baseval); // 0 1 2 3 => G A T C
int base_to_int_value(char nucleotide); // (GATC) = {0 1 2 3}, others = -1

kmer_int_type_t kmer_to_intval(string kmer); // must be less than 32 bases for 64-bit conversion
string decode_kmer_from_intval (kmer_int_type_t intval, unsigned int kmer_length);
kmer_int_type_t revcomp_val(kmer_int_type_t kmer, unsigned int kmer_length);
kmer_int_type_t get_canonical_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length);
kmer_int_type_t get_DS_kmer_val(kmer_int_type_t kmer_val, unsigned int kmer_length);
vector<kmer_int_type_t> sequence_string_to_kmer_int_type_vector(const string& sequence, int kmer_length);



float compute_entropy(string& kmer);
float compute_entropy(kmer_int_type_t kmer, unsigned int kmer_length);

string replace_nonGATC_chars_with_A(string& input_seq);

unsigned long long generateHash(string& s);


#endif


