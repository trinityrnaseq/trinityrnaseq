#include "MaskedKmerCounter.hpp"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <time.h>
#include "stacktrace.hpp"
#include <fstream>
#include <math.h>
#include "sequenceUtil.hpp"
#include "irke_common.hpp"
#include <omp.h>

MaskedKmerCounter::MaskedKmerCounter(unsigned int kmer_length, bool is_ds) {

    if (kmer_length > 32) {
        throw(stacktrace() + "\n\nKmer length exceeds max of 32");  // STORING KMERS AS 64-bit INTEGERS, 2-bit ENCODING.  CANNOT GO HIGHER THAN 32 BASES
    }
    if (kmer_length % 2 != 1) {
        throw(stacktrace() + "\n\nkmer length must be an odd number");
    }
    
    
    this->_kmer_length = kmer_length;
    _DS_MODE = is_ds;
    _masked_pos = kmer_length / 2 + 1;

    _masked_bit_fields = (kmer_int_type_t) 3;
    _masked_bit_fields  = _masked_bit_fields << (_masked_pos * 2 - 2);
    
}


unsigned int MaskedKmerCounter::get_kmer_length() {
    return(_kmer_length);
}

unsigned long MaskedKmerCounter::size() {
    return(_masked_kmer_map.size());
}


void MaskedKmerCounter::describe_kmers() {
    
    Masked_kmer_map_iterator it;
    
    for (it = _masked_kmer_map.begin(); it != _masked_kmer_map.end(); it++) {
        
        kmer_int_type_t kmer_val = it->first;

        cout << describe_kmer(kmer_val) << endl;
    }
    
}

void MaskedKmerCounter::dump_kmers_to_file(string& filename) {
    
    ofstream outfile (filename.c_str());
    
    if (! outfile.is_open()) {
        throw (stacktrace() + " Error, cannot write to file: " + filename);
    }
    
    
    
    Masked_kmer_map_iterator it;
    
    for (it = _masked_kmer_map.begin(); it != _masked_kmer_map.end(); it++) {
        
        kmer_int_type_t kmer_val = (it->second).first;
        unsigned int count = (it->second).second;
        
        string kmer = decode_kmer_from_intval(kmer_val, _kmer_length);
        
        outfile << kmer << " " << count << endl;
        
    }
    
    outfile.close();
    
}


Masked_kmer_map_iterator MaskedKmerCounter::find_kmer(kmer_int_type_t kmer_val) {

	if (_DS_MODE) {
        
        kmer_val = get_DS_kmer_val(kmer_val, _kmer_length);
    }

	return _masked_kmer_map.find(kmer_val);
}


bool MaskedKmerCounter::kmer_exists(string kmer) {
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    return(kmer_exists(kmer_val));
    
    
}

bool MaskedKmerCounter::kmer_exists(kmer_int_type_t kmer_val) {
    
	return(get_kmer_count(kmer_val) > 0);

}


unsigned int MaskedKmerCounter::get_kmer_count(string kmer) {
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    return (get_kmer_count(kmer_val));
    
}


unsigned int MaskedKmerCounter::get_kmer_count(kmer_int_type_t kmer_val) {
        
	Masked_kmer_map_iterator it = find_kmer(kmer_val);
    
    if (it != _masked_kmer_map.end())
   		return( (it->second).second);
    else
   		return 0;
    
}

Kmer_Occurrence_Pair MaskedKmerCounter::get_kmer_occurrence_pair(kmer_int_type_t masked_kmer_val) {

    Masked_kmer_map_iterator it = find_kmer(masked_kmer_val);
    
    if (it == _masked_kmer_map.end()) {
        stringstream ss;
        ss << decode_kmer_from_intval(masked_kmer_val, _kmer_length)  << " not found in _masked_kmer_map";
        throw(stacktrace() + "\n\n" + ss.str());
    }
    
    return(it->second);
}


Kmer_Occurrence_Pair MaskedKmerCounter::get_kmer_occurrence_pair(string kmer) {


    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    return(get_kmer_occurrence_pair(kmer_val));
}


string MaskedKmerCounter::get_kmer_string(kmer_int_type_t kmer_val) {
    
    string kmer = decode_kmer_from_intval(kmer_val, _kmer_length);  
    
    return(kmer);
}


kmer_int_type_t MaskedKmerCounter::get_kmer_intval(string kmer) {
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    return(kmer_val);
}


kmer_int_type_t MaskedKmerCounter::get_masked_kmer_intval(kmer_int_type_t kmer_val) {


    //cerr << "incoming: " << decode_kmer_from_intval(kmer_val, _kmer_length) << endl;

    
    kmer_int_type_t masked_kmer_val = kmer_val | _masked_bit_fields;

    //cerr << "masked:   "   << decode_kmer_from_intval(masked_kmer_val, _kmer_length) << endl << endl;
    
    //kmer_int_type_t masked_kmer_val = kmer_val;
    
    return(masked_kmer_val);
}
                                                          

bool MaskedKmerCounter::set_kmer_occurrence_pair(kmer_int_type_t unmasked_kmer_val, unsigned int count) {

    if (_DS_MODE) {
        unmasked_kmer_val = get_DS_kmer_val(unmasked_kmer_val, _kmer_length);
    }
    
    kmer_int_type_t masked_kmer_val = get_masked_kmer_intval(unmasked_kmer_val);
    
    bool store_kmer_occurrence = false; 

    if (kmer_exists(masked_kmer_val)) {
        Kmer_Occurrence_Pair kmer_occurrence_pair = get_kmer_occurrence_pair(masked_kmer_val);
        kmer_int_type_t stored_kmer_val = kmer_occurrence_pair.first;
        unsigned int stored_count = kmer_occurrence_pair.second;

        if (stored_count < count) {
            store_kmer_occurrence = true;
        }
        else if (stored_count == count && unmasked_kmer_val < stored_kmer_val) {
            // count tie situation, store the lower kmer value to be consistent
            store_kmer_occurrence = true;
        } else {
            // leave false
            ;
        }
    } else {
        // haven't seen this masked kmer yet.  store it.
        store_kmer_occurrence = true;
    }
    
    if (store_kmer_occurrence) {
        #pragma omp critical (HashMap)
        {
            _masked_kmer_map[masked_kmer_val] = Kmer_Occurrence_Pair(unmasked_kmer_val, count);
        }
        return(true);
    }
    
    return(false);
    
}

bool MaskedKmerCounter::set_kmer_occurrence_pair(string kmer, unsigned int count) {
    
    // no storing of kmers containing non-GATC characters
    if (contains_non_gatc(kmer)) { 
        return(false);
    }
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    bool ret = set_kmer_occurrence_pair(kmer_val, count);

    return(ret);
}


string MaskedKmerCounter::describe_kmer(string& kmer) {
    
	kmer_int_type_t kmer_val = get_kmer_intval(kmer);
    return(describe_kmer(kmer_val));
}

string MaskedKmerCounter::describe_kmer(kmer_int_type_t masked_kmer_val) {

    string masked_key_kmer = decode_kmer_from_intval(masked_kmer_val, _kmer_length);

    Kmer_Occurrence_Pair kp = get_kmer_occurrence_pair(masked_kmer_val);
    string stored_unmasked_kmer = decode_kmer_from_intval(kp.first, _kmer_length);
    unsigned int count = kp.second;

    stringstream ss;
    ss << "key: " << masked_key_kmer << ", stored unmasked kmer: " << stored_unmasked_kmer << ", with count: " << count;

    return(ss.str());
}


const Masked_kmer_map& MaskedKmerCounter::get_masked_kmer_map() const {
    
    return (_masked_kmer_map);
    
}

bool MaskedKmerCounter::is_DS_mode() {
    return(_DS_MODE);
}
