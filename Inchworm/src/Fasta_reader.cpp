#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "stacktrace.hpp"
#include <algorithm>
#include <omp.h>

//constructor
Fasta_reader::Fasta_reader (string filename) {
  
    //this->_hasNext = false;
    
    this->end_reading = -1; // turn off
    this->file_byte_pos = 0; // init

    if (filename == "-") {
        filename = "/dev/fd/0"; // read from stdin
    }
    this->_filereader.open(filename.c_str());
    if (! _filereader.is_open()) {
        throw(stacktrace() + "\n\nError, cannot open file " + filename );
    }
        
    this->_init_reader();
   
}

Fasta_reader::Fasta_reader(string filename, long start_reading, long end_reading) {

    this->file_byte_pos = start_reading;
    this->end_reading = end_reading;
    
    this->_filereader.open(filename.c_str());
    if (! _filereader.is_open()) {
        throw(stacktrace() + "\n\nError, cannot open file " + filename );
    }

    if (start_reading > 0) {
      this->_filereader.seekg(start_reading, _filereader.beg);
      this->file_byte_pos = start_reading;
    }
    
    this->_init_reader();
    


}

void Fasta_reader::_init_reader() {

    // primer reader to first fasta header
    getline(this->_filereader, this->_lastline);
    this->file_byte_pos += this->_lastline.length() + 1;

    while ((! this->_filereader.eof()) && this->_lastline[0] != '>') {
        getline(this->_filereader, this->_lastline);
        this->file_byte_pos += this->_lastline.length() + 1;
    }

}


bool Fasta_reader::hasNext() {
    bool ret;
    
    #pragma omp critical (FileReader)
    {
        ret = !(this->_filereader.eof());
        
        // see if we're reading only part of the file.
        if (ret && this->end_reading > 0) {
            if (this->file_byte_pos >= end_reading) { // bad:  this->_filereader.tellg() >= end_reading) {
                // force it to go to the end of the file
                this->_filereader.seekg(0, this->_filereader.end);
                ret = false;
            }
        }
    }
    return ret;
}


Fasta_entry Fasta_reader::getNext() {
    
    string sequence;
    string header;
    bool ret;

    #pragma omp critical (FileReader)
    {
        header = this->_lastline;
        
        ret = !(this->_filereader.eof());
        if (ret == true)
        {
            this->_lastline = "";
            while ((! this->_filereader.eof()) && this->_lastline[0] != '>') {
                getline(this->_filereader, this->_lastline);
                
                this->file_byte_pos += this->_lastline.length() + 1;

                // cerr << "Comparing mine: " << this->file_byte_pos << " to " << this->_filereader.tellg() << endl;
                
                if (this->_lastline[0] != '>') {
                    sequence += this->_lastline;
                }
            }
        }

        // check if only reading section of a file
        if (this->end_reading > 0) {
            if (this->file_byte_pos >= end_reading) { // bad: this->_filereader.tellg() >= end_reading) {
                // force it to go to the end of the file
                this->_filereader.seekg(0, this->_filereader.end);
            }
        }

    }
    
    if (ret == true)
    {
        sequence = remove_whitespace(sequence);
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        Fasta_entry fe(header, sequence);
        return(fe);
    } else {
        Fasta_entry fe("", "");
        return(fe);
    }
}

map<string,string> Fasta_reader::retrieve_all_seqs_hash() {
    
    map<string,string> all_seqs_hash;
    
    while (this->hasNext()) {
        Fasta_entry f = this->getNext();
        string acc = f.get_accession();
        string seq = f.get_sequence();
        
        all_seqs_hash[acc] = seq;
    }
    
    return(all_seqs_hash);
}
