#include "analysis/DeBruijnGraph.h"
#include "analysis/string_util.h"
#include <string>
#include <sstream>
#include <queue>
#include <algorithm>


// class vars
unsigned int DeBruijnKmer::id_counter = 0;

const char DeBruijnKmer::_G_mask = 8;
const char DeBruijnKmer::_A_mask = 4;
const char DeBruijnKmer::_T_mask = 2;
const char DeBruijnKmer::_C_mask = 1;

struct DeBruijnKmerSorter {
    bool operator() (const DeBruijnKmer& a, const DeBruijnKmer& b) {
        return (a.getID() < b.getID());
    }
} kmer_sorter;



DeBruijnKmer::DeBruijnKmer(kmer_int_type_t k) { 

    _kmer = k;
    _prev = 0;
    _next = 0;
    
    id = ++id_counter;
}

DeBruijnKmer::DeBruijnKmer(const DeBruijnKmer& dk) { 
    
    _kmer = dk._kmer;
    _prev = dk._prev;
    _next = dk._next;
    
    id = dk.id;
    
}



unsigned int DeBruijnKmer::getID() const {

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

void DeBruijnKmer::add_next_kmer(kmer_int_type_t k, unsigned int kmer_length) {
    
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

}

void DeBruijnGraph::add_sequence(const string& sequence) {

    vector<kmer_int_type_t> kit_vec = sequence_string_to_kmer_int_type_vector(sequence, _kmer_length);
    
    for (int i = 0; i < (int) kit_vec.size(); i++) {
        kmer_int_type_t k = kit_vec[i];
        DeBruijnKmer& dk = get_kmer_node(k);
        if (i != 0) {
            // set prev kmer
            kmer_int_type_t pk = kit_vec[i-1];
            dk.add_prev_kmer(pk, _kmer_length);
        }
        if (i != (int) kit_vec.size()-1) {
            // set next kmer
            kmer_int_type_t nk = kit_vec[i+1];
            dk.add_next_kmer(nk, _kmer_length);
        }
    }
}


DeBruijnKmer& DeBruijnGraph::get_kmer_node(kmer_int_type_t k) {

    DeBruijnKmerMap::iterator it = _kmer_map.find(k);
    
    if (it == _kmer_map.end()) {
        // add it
        DeBruijnKmer dk(k);
        _kmer_map.insert(pair<kmer_int_type_t,DeBruijnKmer>(k, dk));
        it = _kmer_map.find(k);
    }
    return(it->second);
    
}


string DeBruijnGraph::toString() {

    stringstream s;

    vector<DeBruijnKmer> dk_vec;

    for (DeBruijnKmerMap::iterator it = _kmer_map.begin();
         it != _kmer_map.end();
         it++) {

        DeBruijnKmer& dk = it->second;
        dk_vec.push_back(dk);
    }

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
    
    vector<DeBruijnKmer> dk_vec;

    if (! sStrand) {
        // got DS mirror graph
        dk_vec = deconvolute_DS_mirror_graph();
    }
    else {
        for (DeBruijnKmerMap::iterator it = _kmer_map.begin();
             it != _kmer_map.end();
             it++) {

            dk_vec.push_back(it->second);
        }
    }
    
    // sort vector by id
        
    sort(dk_vec.begin(), dk_vec.end(), kmer_sorter);
    
    // create Chrysalis graph format
    
    for (int i = 0; i < (int) dk_vec.size(); i++) {
        
        DeBruijnKmer& dk = dk_vec[i];
        unsigned int node_id = dk.getID();
        string kmer = decode_kmer_from_intval(dk.get_kmer_int_val(), _kmer_length);

        vector<kmer_int_type_t> prev_kmers = dk.get_prev_kmers(_kmer_length);
        
        if ( ! prev_kmers.empty() ) {
            for (int j = 0; j < (int) prev_kmers.size(); j++) {
                kmer_int_type_t k = prev_kmers[j];
                DeBruijnKmer& pdk = _kmer_map.find(k)->second;
                unsigned int p_id = pdk.getID();
                
                s << node_id << "\t" << p_id << "\t" << 1 << "\t" << kmer << "\t" << 1 << endl;
            }
        }
        else {
            // root kmer
            s << node_id << "\t" << -1 << "\t" << 1 << "\t" << kmer << "\t" << 1 << endl;
        }
    }
    
    
    
    return(s.str());
}

vector<kmer_int_type_t> DeBruijnGraph::get_root_kmers() {

    vector<kmer_int_type_t> kit_vec;

    for (DeBruijnKmerMap::iterator it = _kmer_map.begin();
         it != _kmer_map.end();
         it++) {
        
        DeBruijnKmer& dk = it->second;
        if (dk.get_prev_kmers(_kmer_length).size() == 0) {
            kit_vec.push_back(it->first);
        }
        
    }

    return(kit_vec);
}




vector<DeBruijnKmer> DeBruijnGraph::deconvolute_DS_mirror_graph() {

    // cout << "Deconvoluting DS mirror" << endl;
    
    map<kmer_int_type_t,bool> seen;
    

    vector<DeBruijnKmer> dk_vec;
    
    vector<kmer_int_type_t> root_kmers = get_root_kmers();
    if (root_kmers.size() == 0) {
        // circular, initiate from the first kmer
        // cerr << "Missing root kmers, taking first kmer node instead." << endl;
        root_kmers.push_back(_kmer_map.begin()->first);
    }


    queue<kmer_int_type_t> kmer_queue;

    for (unsigned int r = 0; r < root_kmers.size(); r++) {

        vector<kmer_int_type_t> collected_kmers;

        kmer_int_type_t rk = root_kmers[r];
        DeBruijnKmer rdk = _kmer_map.find(rk)->second;
        
        // cerr << "Got root node ID: " << rdk.getID() << endl;
        
        kmer_queue.push(rk);
        
        while (! kmer_queue.empty()) {

            kmer_int_type_t k = kmer_queue.front();
            kmer_queue.pop();
            
            if (seen.find(k) != seen.end()) {
                continue;
            }
            
            DeBruijnKmer dk = _kmer_map.find(k)->second;
            dk_vec.push_back(dk);
            seen[k] = true;
            collected_kmers.push_back(k);
            // cerr << "-adding: " << dk.getID() << endl;
            
            // populate overlapping kmers into the queue
            vector<kmer_int_type_t> prev_kmers = dk.get_prev_kmers(_kmer_length);
            for (int i = 0; i < (int) prev_kmers.size(); i++) {
                kmer_int_type_t k = prev_kmers[i];
                if (seen.find(k) == seen.end()) {
                    // cerr << "adding prev kmer: " << k << endl;
                    kmer_queue.push(k);
                }
            }
            // now check forward kmers
            vector<kmer_int_type_t> next_kmers = dk.get_next_kmers(_kmer_length);
            for (int i = 0; i < (int) next_kmers.size(); i++) {
                kmer_int_type_t k = next_kmers[i];
                if (seen.find(k) == seen.end()) {
                    // cerr << "adding next kmer: " << k << endl;
                    kmer_queue.push(k);
                }
            }
                        
        }
        
        // set revcomp of collected kmers to found
        for (int i = 0; i < (int) collected_kmers.size(); i++) {
            kmer_int_type_t k = collected_kmers[i];
            k = revcomp_val(k, _kmer_length);
            seen[k] = true;
        }

        

    }
    

    return(dk_vec);
    

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
    
    for (int i = 0; i < left_extensions.size(); i++) {


        for (int j = 0; j < right_extensions.size(); j++) {

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

    for (int i = 0; i < adjacent_kmers.size(); i++) {
        
        kmer_int_type_t adjacent_kmer = adjacent_kmers[i];
        if (kmer_seen.find(adjacent_kmer) != kmer_seen.end()) {
            // already explored.
            continue;
        }
        
        string kmer_string = decode_kmer_from_intval(adjacent_kmer, get_kmer_length());
        
        char char_extend = (direction == 'R') ? kmer_string[kmer_string.length()-1] : kmer_string[0];

        kmer_extension_chars.push_back(char_extend);
        
        if (kmer_extension_chars.size() == flank_extension_length) {
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

    
