#include "DeBruijnGraph.hpp"
#include "string_util.hpp"
#include <string>
#include <sstream>
#include <queue>
#include <algorithm>
#include "irke_common.hpp"
#include "stacktrace.hpp"


const char DeBruijnKmer::_G_mask = 8;
const char DeBruijnKmer::_A_mask = 4;
const char DeBruijnKmer::_T_mask = 2;
const char DeBruijnKmer::_C_mask = 1;

struct DeBruijnKmerSorter {
    bool operator() (const DeBruijnKmer& a, const DeBruijnKmer& b) {
        return (a.getID() < b.getID());
    }
} kmer_sorter;



DeBruijnKmer::DeBruijnKmer(kmer_int_type_t k, long long kmer_id) { 

    _kmer = k;
    _prev = 0;
    _next = 0;
    
    
    id = kmer_id;
    
}

DeBruijnKmer::DeBruijnKmer(const DeBruijnKmer& dk) { 
    
    _kmer = dk._kmer;
    _prev = dk._prev;
    _next = dk._next;
    
    id = dk.id;
    
}


string DeBruijnKmer::toString(int kmer_length) {
    
    stringstream text;

    text << "[" << getID() << "\t" << (int) _prev << "\t" << decode_kmer_from_intval(_kmer, kmer_length) << "\t" << (int) _next << "]";
    
    return(text.str());
}


long long DeBruijnKmer::getID() const {

    return(id);
}

kmer_int_type_t DeBruijnKmer::get_kmer_int_val() const {
    
    return(_kmer);
}



vector<kmer_int_type_t> DeBruijnKmer::get_prev_kmers(unsigned int kmer_length) {

    vector<kmer_int_type_t> prev_kmers;
    
    kmer_int_type_t reverse_suffix = _kmer >> 2;
    
    // examine G
    if (_prev & _G_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('G');
        k = (k << (kmer_length*2-2)) | reverse_suffix;
        prev_kmers.push_back(k);
    }
    
    // examine A
    if (_prev & _A_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('A');
        k = (k << (kmer_length*2-2)) | reverse_suffix;
        prev_kmers.push_back(k);
    }
    
    // examine T
    if (_prev & _T_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('T');
        k = (k << (kmer_length*2-2)) | reverse_suffix;
        prev_kmers.push_back(k);
    }

    // examine C
    if (_prev & _C_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('C');
        k = (k << (kmer_length*2-2)) | reverse_suffix;
        prev_kmers.push_back(k);
    }

    return(prev_kmers);
}



void DeBruijnKmer::add_prev_kmer(kmer_int_type_t k, unsigned int kmer_length) {


    if (IRKE_COMMON::MONITOR >= 3) {
        string prev_kmer = decode_kmer_from_intval(k, kmer_length);
        cerr << "Adding prev_kmer: " << prev_kmer << " to kmer: " << decode_kmer_from_intval(_kmer, kmer_length) << endl;
    }
    

    k = (k >> (kmer_length*2-2)); // get the first base of the previous kmer


    
    if (k == (kmer_int_type_t) base_to_int_value('G')) {
        _prev |= _G_mask;
    }
    if (k == (kmer_int_type_t) base_to_int_value('A')) {
        _prev |= _A_mask;
    }
    if (k == (kmer_int_type_t) base_to_int_value('T')) {
        _prev |= _T_mask;
    }
    if (k == (kmer_int_type_t) base_to_int_value('C')) {
        _prev |= _C_mask;
    }

}


vector<kmer_int_type_t> DeBruijnKmer::get_next_kmers(unsigned int kmer_length) {

    vector<kmer_int_type_t> next_kmers;
    
    kmer_int_type_t forward_prefix = (_kmer << (33-kmer_length)*2) >> (32-kmer_length)*2;

    // examine G
    if (_next & _G_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('G');
        k |= forward_prefix;
        next_kmers.push_back(k);
    }
    
    // examine A
    if (_next & _A_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('A');
        k |= forward_prefix;
        next_kmers.push_back(k);
    }
    
    // examine T
    if (_next & _T_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('T');
        k |= forward_prefix;
        next_kmers.push_back(k);
    }

    // examine C
    if (_next & _C_mask) {
        kmer_int_type_t k = (kmer_int_type_t) base_to_int_value('C');
        k |= forward_prefix;
        next_kmers.push_back(k);
    }

    return(next_kmers);
}

void DeBruijnKmer::add_next_kmer(kmer_int_type_t k) {
    
    //cerr << "input: " << k;
    
    k &= 3ull; // get the last base of the next kmer
    
    // cerr << ", now: " << k << endl;
    
    if (k == (kmer_int_type_t) base_to_int_value('G')) {
        _next |= _G_mask;
    }
    if (k == (kmer_int_type_t) base_to_int_value('A')) {
        _next |= _A_mask;
    }
    if (k == (kmer_int_type_t) base_to_int_value('T')) {
        _next |= _T_mask;
    }
    if (k == (kmer_int_type_t) base_to_int_value('C')) {
        _next |= _C_mask;
    }

}


DeBruijnGraph::DeBruijnGraph(unsigned int kmer_length) {

    _kmer_length = kmer_length;
    _kmer_id_counter = 0;

}

void DeBruijnGraph::add_sequence(const string& sequence, bool sStrand) {

    vector<kmer_int_type_t> kit_vec = sequence_string_to_kmer_int_type_vector(sequence, _kmer_length);
    
    for (int i = 0; i < (int) kit_vec.size(); i++) {
        kmer_int_type_t k = kit_vec[i];

        kmer_int_type_t orig_sequence_k = k;
        if (! sStrand) {
            k = get_canonical_kmer_val(k, _kmer_length);
        }
        
        DeBruijnKmer& dk = get_kmer_node(k);
        if (i != 0) {
            // set prev kmer
            kmer_int_type_t pk = kit_vec[i-1];

            if (sStrand) {
                // easy, use as is
                dk.add_prev_kmer(pk, _kmer_length);
            }
            else {
                // harder... see which orientation we're in.
                if (k == orig_sequence_k) {
                    // fine, in original orientation.
                    dk.add_prev_kmer(pk, _kmer_length);
                }
                else {
                    // must work in revcomp mode.
                    pk = revcomp_val(pk, _kmer_length);
                    dk.add_next_kmer(pk);
                }
            }
        }
        if (i != (int) kit_vec.size()-1) {
            // set next kmer
            kmer_int_type_t nk = kit_vec[i+1];

            if (sStrand) {
                // easy, use as is
                dk.add_next_kmer(nk);
            }
            else {
                // harder.... see which orientation we're in.
                if (k == orig_sequence_k) {
                    // fine, in original orietation
                    dk.add_next_kmer(nk);
                }
                else {
                    // must work in revcomp model
                    nk = revcomp_val(nk, _kmer_length);
                    dk.add_prev_kmer(nk, _kmer_length);
                }
            }
        }
    }
}


DeBruijnKmer& DeBruijnGraph::get_kmer_node(kmer_int_type_t k) {
    
    #pragma omp critical
    {
        DeBruijnKmerMap::iterator it = _kmer_map.find(k);
        
        if (it == _kmer_map.end()) {
            // add it
            DeBruijnKmer dk(k, ++_kmer_id_counter);
            _kmer_map.insert(pair<kmer_int_type_t,DeBruijnKmer>(k, dk));
          
            if (IRKE_COMMON::MONITOR >= 3) 
                cerr << "Adding Node: " << dk.toString(get_kmer_length()) << " to graph." << ", now kmer_map is size: " << _kmer_map.size() << endl;
            
            
        }

    }

    DeBruijnKmerMap::iterator it = _kmer_map.find(k);
    return(it->second);
    

}


string DeBruijnGraph::toString() {

    stringstream s;

    vector<DeBruijnKmer> dk_vec;

    for (DeBruijnKmerMap::iterator it = _kmer_map.begin();
         it != _kmer_map.end();
         it++) {

        DeBruijnKmer dk = it->second;
        dk_vec.push_back(dk);
    }
    
    
    if (IRKE_COMMON::MONITOR >= 2) 
        cerr << "Kmer map has size: " << _kmer_map.size() << ", and vector has size: " << dk_vec.size() << endl;
    
    sort(dk_vec.begin(), dk_vec.end(), kmer_sorter);

    for (int i = 0; i < (int) dk_vec.size(); i++) {
        
        DeBruijnKmer& dk = dk_vec[i];

        kmer_int_type_t k = dk.get_kmer_int_val();
        
        string kmer = decode_kmer_from_intval(k, _kmer_length);
        s << "node:" << dk.getID() << "\t" 
          << (int) dk._prev << "\t" 
          << kmer << "\t" 
          << (int) dk._next << endl;
    }
    
    return(s.str());
}




string DeBruijnGraph::toChrysalisFormat(int component_id, bool sStrand) {

    stringstream s;

    s << "Component " << component_id << endl;
    
    //vector<DeBruijnKmer> dk_vec;
    map<kmer_int_type_t,bool> seen;


    vector<kmer_int_type_t> root_kmers = get_root_kmers(sStrand);
    if (root_kmers.size() == 0) {
        // circular, initiate from the first kmer
        // cerr << "Missing root kmers, taking first kmer node instead." << endl;
        root_kmers.push_back(_kmer_map.begin()->first);
    }


    queue<kmer_int_type_t> kmer_queue; //kmers must be oriented to '+' here.
    
    for (unsigned int r = 0; r < root_kmers.size(); r++) {
            
        kmer_int_type_t rk = root_kmers[r];  // should always be considered as the starting '+' orientation
        
        // cerr << "Got root node ID: " << rdk.getID() << endl;
        
        kmer_queue.push(rk); 
        
        while (! kmer_queue.empty()) {

            kmer_int_type_t k = kmer_queue.front();  // in the relevant orientation.
            kmer_queue.pop();
                        
            if (seen.find(k) != seen.end()) {
                // already seen it.
                continue;
            }
            seen[k] = true;
            if (! sStrand) {
                seen[ revcomp_val(k, _kmer_length) ] = true;
            }
            
            string kmer_seq = decode_kmer_from_intval(k, _kmer_length);
            
            
            // get the dk object:
            // will depend on sStrand or not as to what the kmer represnetation is there.
            kmer_int_type_t dk_stored_kmer_val = (sStrand) ? k : get_canonical_kmer_val(k, _kmer_length);

            char dk_orient = (dk_stored_kmer_val == k) ? '+' : '-';
            
            DeBruijnKmer dk = _kmer_map.find(dk_stored_kmer_val)->second;
            long long dk_id = dk.getID();

            /////////////////////////////////////////////
            // populate overlapping kmers into the queue

            /////////////////////////
            // walk kmer to the left:
            
            vector<kmer_int_type_t> prev_kmers;
            if (dk_orient == '+') {
                prev_kmers = dk.get_prev_kmers(_kmer_length);
            }
            else {
                // must get the next kmers, and revcomp them so they're in the correct orientation.
                vector<kmer_int_type_t> tmp_prev_kmers = dk.get_next_kmers(_kmer_length);
                for (int i=0; i < (int) tmp_prev_kmers.size(); i++) {
                    prev_kmers.push_back(revcomp_val( tmp_prev_kmers[i], _kmer_length));
                }
            }


            // kmer graph reporting.
            if (prev_kmers.size() > 0) {
                for (int i=0; i < (int) prev_kmers.size(); i++) {
                    
                    kmer_int_type_t pk_kmer = prev_kmers[i];
                    kmer_int_type_t pk_stored_kmer = (sStrand) ? pk_kmer : get_canonical_kmer_val(pk_kmer, _kmer_length);
                    DeBruijnKmer pkd = _kmer_map.find(pk_stored_kmer)->second; 
                    long long pkd_id = pkd.getID();
                    //string pk_kmer_seq = decode_kmer_from_intval(pk_kmer, _kmer_length);
                    
                    s << dk_id << "\t" << pkd_id << "\t" << 1 << "\t" << kmer_seq << "\t" << 1 << endl;
                    
                    kmer_queue.push(pk_kmer);
                }
            }
            else {
                // hit left end, no prev extension.
                s << dk_id << "\t" << -1 << "\t" << 1 << "\t" << kmer_seq << "\t" << 1 << endl;
            }


            /////////////////////////
            // walk kmer to the right:

            // just add any next kmers to the queue.
            
            vector<kmer_int_type_t> next_kmers;
            if (dk_orient == '+') {
                next_kmers = dk.get_next_kmers(_kmer_length);
            }
            else {
                // must get the prev kmers, and revcomp them so they're in the correct orientation.
                vector<kmer_int_type_t> tmp_next_kmers = dk.get_prev_kmers(_kmer_length);
                for (int i=0; i < (int) tmp_next_kmers.size(); i++) {
                    next_kmers.push_back(revcomp_val( tmp_next_kmers[i], _kmer_length));
                }
            }
            
            if (next_kmers.size() > 0) {
                for (int i=0; i < (int) next_kmers.size(); i++) {
                    
                    kmer_int_type_t nk_kmer = next_kmers[i];
                    //kmer_int_type_t nk_stored_kmer = (sStrand) ? nk_kmer : get_canonical_kmer_val(nk_kmer, _kmer_length);
                    //DeBruijnKmer nkd = _kmer_map.find(nk_stored_kmer)->second; 
                    //long long nkd_id = nkd.getID();
                    //string nk_kmer_seq = decode_kmer_from_intval(nk_kmer, _kmer_length);
                    
                    //s << nkd_id << "\t" << dk_id << "\t" << 1 << "\t" << nk_kmer_seq << "\t" << 1 << endl;
                    
                    kmer_queue.push(nk_kmer);
                }
            }
                        
        } // end of kmer queue

    } // end of root kmer traversal

    
    return(s.str());
}

vector<kmer_int_type_t> DeBruijnGraph::get_root_kmers(bool sStrand) {

    // start nodes
    // or end nodes in the graph if not strand-specific
    
    vector<kmer_int_type_t> kit_vec;

    for (DeBruijnKmerMap::iterator it = _kmer_map.begin();
         it != _kmer_map.end();
         it++) {
        
        DeBruijnKmer& dk = it->second;
        if (dk.get_prev_kmers(_kmer_length).size() == 0
            ||
            ( (!sStrand) && dk.get_next_kmers(_kmer_length).size() == 0) // ds could go either way
            ) {
            kit_vec.push_back(it->first);
        }
        
    }

    return(kit_vec);
}


bool DeBruijnGraph::kmerExists(kmer_int_type_t kval) {

    if (_kmer_map.find(kval) == _kmer_map.end()) {
        return(false);
    }
    else {
        return(true);
    }
}

unsigned int DeBruijnGraph::get_kmer_length() {
    return(_kmer_length);
}



string DeBruijnGraph::get_kmer_from_char_vector(vector<char> char_vec, char direction) {

    if (direction == 'L') {
        reverse(char_vec.begin(), char_vec.end());
    }

    stringstream kmer_string;
    for (int i = 0; i < (int) char_vec.size(); i++) {
        kmer_string << char_vec[i];
    }

    return(kmer_string.str());
}


vector<string> DeBruijnGraph::get_candidate_weldmers(kmer_int_type_t kmer_val, int weldmer_length) {
    
    string kmer_string = decode_kmer_from_intval(kmer_val, get_kmer_length());
    
    
    // cerr << "Graph seed for weldmer construction: " << kmer_string << endl;

    vector<string> right_extensions;
    vector<string> left_extensions;
    
    //----------------------------------------------------------------------------------------
    //-----  Create combinations of <left--kmer--right> sequences as candidate weldmers ------
    //----------------------------------------------------------------------------------------
    
    int flank_extension_length = (weldmer_length - get_kmer_length())/2;
    
    vector<char> extension_vector;
    map<kmer_int_type_t,bool> kmer_seen;

    recursively_construct_kmer_extensions(kmer_val, extension_vector, right_extensions, 'R', kmer_seen, flank_extension_length);

    recursively_construct_kmer_extensions(kmer_val, extension_vector, left_extensions, 'L', kmer_seen, flank_extension_length);

    // join them up:
    
    vector<string> weldmers;
    
    for (size_t i = 0; i < left_extensions.size(); i++) {


        for (size_t j = 0; j < right_extensions.size(); j++) {

            string weldmer = left_extensions[i] + kmer_string + right_extensions[j];

            weldmers.push_back(weldmer);
        }
    }
    
    return(weldmers);
    
}


void DeBruijnGraph::recursively_construct_kmer_extensions(kmer_int_type_t kmer_val,
                                                          vector<char>& kmer_extension_chars,
                                                          vector<string>& extension_kmer_strings,
                                                          char direction, // L | R  fg left or right from the seed.
                                                          map<kmer_int_type_t,bool>& kmer_seen,
                                                          int flank_extension_length
                                                          )  {
    
    kmer_seen[kmer_val] = true;
    
    DeBruijnKmer& kmer_node = get_kmer_node(kmer_val);
    
    // cerr << kmer_extension_chars.size() << " vs " << flank_extension_length << ", depth, " << direction << ", " << decode_kmer_from_intval(kmer_val, get_kmer_length()) << endl;
    
    vector<kmer_int_type_t> adjacent_kmers;
    
    if (direction == 'R') { // go right
        adjacent_kmers = kmer_node.get_next_kmers(get_kmer_length());
    }
    else {
        // go left
        adjacent_kmers = kmer_node.get_prev_kmers(get_kmer_length());
    }

    for (size_t i = 0; i < adjacent_kmers.size(); i++) {
        
        kmer_int_type_t adjacent_kmer = adjacent_kmers[i];
        if (kmer_seen.find(adjacent_kmer) != kmer_seen.end()) {
            // already explored.
            continue;
        }
        
        string kmer_string = decode_kmer_from_intval(adjacent_kmer, get_kmer_length());
        
        char char_extend = (direction == 'R') ? kmer_string[kmer_string.length()-1] : kmer_string[0];

        kmer_extension_chars.push_back(char_extend);
        
        if (kmer_extension_chars.size() == (size_t)flank_extension_length) {
            // reached maximal extension.
            string kmer_extension_string = get_kmer_from_char_vector(kmer_extension_chars, direction);
            extension_kmer_strings.push_back(kmer_extension_string);
        }
        else {
            // continue the recursive extension
            recursively_construct_kmer_extensions(adjacent_kmer, kmer_extension_chars, extension_kmer_strings, direction, kmer_seen, flank_extension_length);
        }
        
        kmer_extension_chars.pop_back();
        kmer_seen.erase(adjacent_kmer);
    }
    
    return;
        
    
}

    
