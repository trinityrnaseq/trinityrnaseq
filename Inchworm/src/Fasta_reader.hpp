#ifndef __FASTA_READER__
#define __FASTA_READER__

#include <iostream>
#include <fstream>
#include <string>
#include "Fasta_entry.hpp"
#include <map>

using namespace std;

class Fasta_reader {
    
public:
    Fasta_reader(string filename); // constructor
    Fasta_reader(string filename, long start_reading, long end_reading); // constructor
    
    bool hasNext();  // prime for next 
    Fasta_entry getNext(); // retrieves next fasta entry, or NULL if EOF.
    
    map<string,string> retrieve_all_seqs_hash();
    
    
    
private:
    ifstream _filereader;
    //bool _hasNext;
    string _lastline;
    long end_reading; // optional file position to stop reading.

    long file_byte_pos;


    void _init_reader();

};


#endif
