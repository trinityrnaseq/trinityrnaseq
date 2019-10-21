#include "KmerCounter.hpp"
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

KmerCounter::KmerCounter(unsigned int kmer_length, bool is_ds) {

    if (kmer_length > 32) {
        throw(stacktrace() + "\n\nKmer length exceeds max of 32");  // STORING KMERS AS 64-bit INTEGERS, 2-bit ENCODING.  CANNOT GO HIGHER THAN 32 BASES
    }
    
    this->_kmer_length = kmer_length;
    _DS_MODE = is_ds;
    
    
#ifdef __GOOGLE__
    
    kmer_int_type_t delete_val = 0xFFFFFFFFFFFFFFFFll;
    
    _kmer_counter.set_deleted_key(delete_val);
    
    
#endif
    
}

void KmerCounter::add_sequence(string& sequence, unsigned int cov) {
    
    for (unsigned int i = 0; i <= sequence.length() - _kmer_length; i++) {
        
        string kmer = sequence.substr(i, _kmer_length);
        
        add_kmer(kmer, cov);  // if DS-mode, kmer is added in both orientations

    }
	
}


unsigned int KmerCounter::get_kmer_length() {
    return(_kmer_length);
}

unsigned long KmerCounter::size() {
    return(_kmer_counter.size());
}


void KmerCounter::describe_kmers() {
    
    Kmer_counter_map_iterator it;
    
    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        
        kmer_int_type_t kmer_val = it->first;
        
        string kmer = decode_kmer_from_intval(kmer_val, _kmer_length);
        
        cout << describe_kmer(kmer) << endl;
        
    }
    
}


bool KmerCounter::prune_kmer_extensions( float min_ratio_non_error) {
    
    Kmer_counter_map_iterator it;
    vector<kmer_int_type_t> deletion_list;
    
    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        
        kmer_int_type_t kmer_val = it->first;
        int count = it->second;
        if (count == 0) {
            continue;
        }
        
        string kmer_str = get_kmer_string(kmer_val); 
        

    	vector<Kmer_Occurence_Pair> candidates = get_forward_kmer_candidates(kmer_val);
        
        if (candidates.size() > 1) {
            
            int dominant_count = candidates[0].second;
            
            for (unsigned int i = 1; i < candidates.size(); i++)
                {
                    if (candidates[i].second)
                        {
                            int candidate_count = candidates[i].second;
                            float ratio_dominant_count = (float) candidate_count/dominant_count;
                            
                            if (dominant_count > 0 && ratio_dominant_count < min_ratio_non_error) {
                                
                                Kmer_counter_map_iterator kmer_candidate = find_kmer(candidates[i].first);
                                
                                if (IRKE_COMMON::MONITOR >= 2) {
                                    string kmer_ext_str = get_kmer_string(kmer_candidate->first);
                                    cerr << "Pruning kmer: " << kmer_ext_str << " extension of: " <<  describe_kmer(kmer_str)
                                         << " with ratio dominant count: " << ratio_dominant_count << endl;
                    }
                                
                                deletion_list.push_back(kmer_candidate->first);
                                kmer_candidate->second = 0; // disable when encountered in further iterations.
                                
                            }
                        }
                }
        }
    }
    
    if (deletion_list.size() > 0) {
        
        for (unsigned int i=0; i<deletion_list.size(); i++) {
            prune_kmer(deletion_list[i]);
        }
        
        return(true);
    }
    else {
        return(false);
    }
}


bool KmerCounter::prune_some_kmers(unsigned int min_count, float min_entropy, bool prune_error_kmers, float min_ratio_non_error) {

    Kmer_counter_map_iterator it;
    vector<kmer_int_type_t> deletion_list;

    int count_pruned = 0;

        
    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++)
    {
        kmer_int_type_t kmer = it->first;
        unsigned int count = it->second;

        if (count == 0)
        	continue;

        string kmer_str = get_kmer_string(kmer);

        if (count < min_count)
        {
        	
        	it->second = 0; // immediate removal
        	count_pruned++;
            continue;
        }

        if (compute_entropy(kmer,_kmer_length) < min_entropy)
        {
        	
        	it->second = 0;  // immediate removal
        	count_pruned++;
            continue;
        }
        
        if (prune_error_kmers)
        {

            // examine forward extension kmers
            vector<Kmer_Occurence_Pair> candidates = get_forward_kmer_candidates(kmer); //sorted descendingly
            
            if (candidates.size() > 1) {
                
                int dominant_count = candidates[0].second;
        
                
        
                for (unsigned int i = 1; i < candidates.size(); i++)
                    {
                        if (candidates[i].second)
                            {
                                int candidate_count = candidates[i].second;
                                
                                float ratio_dominant_count = (float) candidate_count/dominant_count;
                                
                                float ratio_curr_count = (float) candidate_count / count;
     
                                if (dominant_count > 0 
                                    && 
                                    ratio_dominant_count < min_ratio_non_error
                                    &&
                                    ratio_curr_count < min_ratio_non_error 
                                    ) {
                                    
                                    Kmer_counter_map_iterator kmer_candidate = find_kmer(candidates[i].first);
                                    
                                    if (IRKE_COMMON::MONITOR >= 2) {
                                        string kmer_ext_str = get_kmer_string(kmer_candidate->first);
                                        cerr << "Pruning kmer: " << kmer_ext_str << " FORWARD extension of: " <<  describe_kmer(kmer_str)
                                             << " with ratio dominant count: " << ratio_dominant_count << endl;
                                    }
                                    
                                    deletion_list.push_back(kmer_candidate->first);
                                    //kmer_candidate->second = 0; // disable when encountered in further iterations.
                                    count_pruned++;
                                    
                                    
                                }
                            }
                    }
            } // end forward kmer pruning.



            // examine reverse extension kmers
            candidates = get_reverse_kmer_candidates(kmer); //sorted descendingly
            
            if (candidates.size() > 1) {
                
                int dominant_count = candidates[0].second;
                
                for (unsigned int i = 1; i < candidates.size(); i++)
                    {
                        if (candidates[i].second)
                            {
                                int candidate_count = candidates[i].second;
                                
                                float ratio_dominant_count = (float) candidate_count/dominant_count;

                                float ratio_curr_count = (float) candidate_count / count;
                                
                                if (dominant_count > 0 
                                    && 
                                    ratio_dominant_count < min_ratio_non_error
                                    &&
                                    ratio_curr_count < min_ratio_non_error
                                    ) {
                                    
                                    Kmer_counter_map_iterator kmer_candidate = find_kmer(candidates[i].first);
                                    
                                    if (IRKE_COMMON::MONITOR >= 2) {
                                        string kmer_ext_str = get_kmer_string(kmer_candidate->first);
                                        cerr << "Pruning kmer: " << kmer_ext_str << " REVERSE extension of: " <<  describe_kmer(kmer_str)
                                             << " with ratio dominant count: " << ratio_dominant_count << endl;
                                    }
                                    
                                    deletion_list.push_back(kmer_candidate->first);
                                    //kmer_candidate->second = 0; // disable when encountered in further iterations.
                                    count_pruned++;
                                    
                                    
                                }
                            }
                    }
            } // end reverse kmer pruning.
            
                        
        }
    }
    
    cerr << "Pruned " << count_pruned << " kmers from catalog." << endl;
    
        
    if (count_pruned > 0) { 
        
        for (unsigned int i=0; i<deletion_list.size(); i++) {
            //prune_kmer(deletion_list[i]);  // dont bother shrinking hashtable
            Kmer_counter_map_iterator kmer_candidate = find_kmer(deletion_list[i]); 
            kmer_candidate->second = 0; // disable when encountered in further iterations.
        }
        
        return(true);
    }
    else {
        return(false);
    }
}


bool KmerCounter::prune_branched_kmers() {

    Kmer_counter_map_iterator it;
    vector<kmer_int_type_t> deletion_list;

    int count_pruned = 0;

    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++)
    {
        kmer_int_type_t kmer = it->first;
        unsigned int count = it->second;

        if (count == 0)
        	continue;

        // examine forward extension kmers
        vector<Kmer_Occurence_Pair> forward_candidates = get_forward_kmer_candidates(kmer); //sorted descendingly
        vector<Kmer_Occurence_Pair> reverse_candidates = get_reverse_kmer_candidates(kmer); //sorted descendingly
        
        if (forward_candidates.size() > 1 || reverse_candidates.size() > 1) {
            // branched kmer
            
            if (IRKE_COMMON::MONITOR >= 2) {
                string kmer_str = get_kmer_string(kmer);
                
                cerr << "Pruning kmer: " << kmer_str << " as branched kmer " << endl;
            }

            //it->second = 0;
            deletion_list.push_back(kmer);
            count_pruned++;
                        
        }

    }
    
    cerr << "Pruned " << count_pruned << " branched kmers from catalog." << endl;
    
    if (count_pruned > 0) {  //deletion_list.size() > 0) {
        
        for (unsigned int i=0; i<deletion_list.size(); i++) {
            //prune_kmer(deletion_list[i]);  // dont bother shrinking hashtable
            Kmer_counter_map_iterator kmer_candidate = find_kmer(deletion_list[i]); 
            kmer_candidate->second = 0; // disable when encountered in further iterations.
        }   
        
        return(true);
    }
    else {
        return(false);
    }
}




void KmerCounter::dump_kmers_to_file(string& filename) {
    
    ofstream outfile (filename.c_str());
    
    if (! outfile.is_open()) {
        throw (stacktrace() + " Error, cannot write to file: " + filename);
    }
    
    
    
    Kmer_counter_map_iterator it;
    
    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        
        kmer_int_type_t kmer_val = it->first;
        unsigned int count = it->second;
        
        string kmer = decode_kmer_from_intval(kmer_val, _kmer_length);
        
        outfile << kmer << " " << count << endl;
        
    }
    
    outfile.close();
    
}


Kmer_counter_map_iterator KmerCounter::find_kmer(kmer_int_type_t kmer_val) {

	if (_DS_MODE) {
        
        kmer_val = get_DS_kmer_val(kmer_val, _kmer_length);
    }

	return _kmer_counter.find(kmer_val);
}


bool KmerCounter::kmer_exists(string kmer) {
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    return(kmer_exists(kmer_val));
    
    
}

bool KmerCounter::kmer_exists(kmer_int_type_t kmer_val) {
    
	return(get_kmer_count(kmer_val) > 0);

}


bool KmerCounter::prune_kmer(string kmer) {
    
        kmer_int_type_t kmer_val = kmer_to_intval(kmer);
        return(prune_kmer(kmer_val));
    
}


bool KmerCounter::prune_kmer(kmer_int_type_t kmer_val) {
    
	bool ret = false;

    #pragma omp critical (HashMap)
    {
      Kmer_counter_map_iterator it = find_kmer(kmer_val);
      if (it != _kmer_counter.end())
	{
	  _kmer_counter.erase(it);
	  ret = true;
	}
    }

    return ret;
}


bool KmerCounter::clear_kmer(kmer_int_type_t kmer_val) {

    bool ret = false;

    // removing the omp critical here to improve performance. thread protection shouldn't be overly important here (bhaas, Jan-5-2014)
    //#pragma omp critical (HashMap)
    {
		Kmer_counter_map_iterator it = find_kmer(kmer_val);
    	if (it != _kmer_counter.end())
    	{
    		it->second = 0;
    		ret = true;
    	}
    }
    
    return ret;
}


unsigned int KmerCounter::get_kmer_count(string kmer) {
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    return (get_kmer_count(kmer_val));
    
}


unsigned int KmerCounter::get_kmer_count(kmer_int_type_t kmer_val) {
        
	Kmer_counter_map_iterator it = find_kmer(kmer_val);
    
    if (it != _kmer_counter.end())
   		return(it->second);
    else
   		return 0;
    
}


string KmerCounter::get_kmer_string(kmer_int_type_t kmer_val) {
    
    string kmer = decode_kmer_from_intval(kmer_val, _kmer_length);  
    
    return(kmer);
}


kmer_int_type_t KmerCounter::get_kmer_intval(string kmer) {
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    return(kmer_val);
}


bool KmerCounter::add_kmer (kmer_int_type_t kmer_val, unsigned int count) {

    if (_DS_MODE) {
        kmer_val = get_DS_kmer_val(kmer_val, _kmer_length);
    }
    
    #pragma omp critical (HashMap)
    {
        _kmer_counter[kmer_val]+= count;
    }
    
    return(true);
    
}



bool KmerCounter::add_kmer(string kmer, unsigned int count) {
    
    // no storing of kmers containing non-GATC characters
    if (contains_non_gatc(kmer)) { 
        return(false);
    }
    
    kmer_int_type_t kmer_val = kmer_to_intval(kmer);
    
    bool ret = add_kmer(kmer_val, count);

    return(ret);
}

vector<Kmer_Occurence_Pair> KmerCounter::get_forward_kmer_candidates(kmer_int_type_t seed_kmer) {

    vector<Kmer_Occurence_Pair> candidates = get_forward_kmer_candidates_unsorted(seed_kmer, false);

    Sort_kmer_by_count_desc sorter (this);
    sort (candidates.begin(), candidates.end(), sorter);

    return(candidates);
}


vector<Kmer_Occurence_Pair> KmerCounter::get_forward_kmer_candidates_unsorted(kmer_int_type_t seed_kmer, bool getZeros) {

	kmer_int_type_t forward_prefix = (seed_kmer << (33-_kmer_length)*2) >> (32-_kmer_length)*2;

    vector<Kmer_Occurence_Pair> candidates;
	for (kmer_int_type_t i=0; i<4; i++)
	{
		Kmer_Occurence_Pair candidate;
		candidate.first = forward_prefix | i;
		candidate.second = get_kmer_count(candidate.first);
		if (candidate.second || getZeros)
			candidates.push_back(candidate);
	}

    return(candidates);
}


kmer_int_type_t*  KmerCounter::get_forward_kmer_candidates_noLookup(kmer_int_type_t seed_kmer, kmer_int_type_t forward_candidates []) {

	kmer_int_type_t forward_prefix = (seed_kmer << (33-_kmer_length)*2) >> (32-_kmer_length)*2;

	for (kmer_int_type_t i=0; i<4; i++) {
		Kmer_Occurence_Pair candidate;
		forward_candidates[i] = forward_prefix | i;
	}

    return(forward_candidates);
}


vector<Kmer_Occurence_Pair> KmerCounter::get_reverse_kmer_candidates(kmer_int_type_t seed_kmer) {

    vector<Kmer_Occurence_Pair> candidates = get_reverse_kmer_candidates_unsorted(seed_kmer, false);

    Sort_kmer_by_count_desc sorter (this);
    sort (candidates.begin(), candidates.end(), sorter);

    return(candidates);
}


vector<Kmer_Occurence_Pair> KmerCounter::get_reverse_kmer_candidates_unsorted(kmer_int_type_t seed_kmer, bool getZeros) {

	kmer_int_type_t reverse_suffix = seed_kmer >> 2;

    vector<Kmer_Occurence_Pair> candidates;
	for (kmer_int_type_t i=0; i<4; i++)
	{
		Kmer_Occurence_Pair candidate;
		candidate.first = (i << (_kmer_length*2-2)) | reverse_suffix;
		candidate.second = get_kmer_count(candidate.first);
		if (candidate.second || getZeros)
			candidates.push_back(candidate);
	}

    return(candidates);
}

kmer_int_type_t* KmerCounter::get_reverse_kmer_candidates_noLookup(kmer_int_type_t seed_kmer, kmer_int_type_t reverse_candidates []) {

	kmer_int_type_t reverse_suffix = seed_kmer >> 2;

	for (kmer_int_type_t i=0; i<4; i++) {
        reverse_candidates[i] = (i << (_kmer_length*2-2)) | reverse_suffix;
	}

    return(reverse_candidates);
}



string KmerCounter::describe_kmer(string& kmer) {
    
	kmer_int_type_t kmer_val = get_kmer_intval(kmer);
	int kmer_count = get_kmer_count(kmer_val);
    
    
    // get counts for reverse candidates:
    vector<Kmer_Occurence_Pair> reverse_candidates = get_reverse_kmer_candidates_unsorted(kmer_val, true);

    int rG_count = reverse_candidates[0].second;
    int rA_count = reverse_candidates[1].second;
    int rT_count = reverse_candidates[2].second;
    int rC_count = reverse_candidates[3].second;
    
    
    // get counts for forward candidates:
    
    vector<Kmer_Occurence_Pair> forward_candidates = get_forward_kmer_candidates_unsorted(kmer_val, true);

    int fG_count = forward_candidates[0].second;
    int fA_count = forward_candidates[1].second;
    int fT_count = forward_candidates[2].second;
    int fC_count = forward_candidates[3].second;
    
    
    stringstream descr;
    
    descr << "G:" << rG_count << " " 
          << "A:" << rA_count << " "
          << "T:" << rT_count << " " 
          << "C:" << rC_count << " "
		
          << "\t" << kmer << ":" << kmer_count << "\t"
		
          << "G:" << fG_count << " "
          << "A:" << fA_count << " "
          << "T:" << fT_count << " " 
          << "C:" << fC_count;
    
    return(descr.str());
    
    
}


void KmerCounter::prune_kmers_min_count(unsigned int min_count) {
    
    Kmer_counter_map_iterator it;
    
    vector<kmer_int_type_t> kmers_to_purge;
    
    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        
        kmer_int_type_t kmer_val = it->first;
        
        unsigned int count = it->second;
        
        if (count < min_count) {
            kmers_to_purge.push_back(kmer_val);
        }
        
    }
    
    for (unsigned int i = 0; i < kmers_to_purge.size(); i++) {
        
        prune_kmer( kmers_to_purge[i] );
    }
    
    return;
}



void KmerCounter::prune_kmers_min_entropy(float min_entropy) {
    
    Kmer_counter_map_iterator it;
    
    vector<kmer_int_type_t> kmers_to_purge;
    
    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {
        
        kmer_int_type_t kmer_val = it->first;
        
        float entropy = compute_entropy(kmer_val,_kmer_length);
        
        if (entropy < min_entropy) {
            kmers_to_purge.push_back(kmer_val);
        }
        
    }
    
    for (unsigned int i = 0; i < kmers_to_purge.size(); i++) {
        
        prune_kmer( kmers_to_purge[i] );
    }
    
    return;
}







const Kmer_counter_map& KmerCounter::get_kmer_counter_map() const {
    
    return (_kmer_counter);
    
}

vector<Kmer_Occurence_Pair> KmerCounter::get_kmers_sort_descending_counts() {
    
    // cerr << "Getting vec of kmers" << endl;
    
    unsigned long start, end;
    unsigned long num_kmers = _kmer_counter.size();

    start = time(NULL);

    //vector<Kmer_counter_map_iterator> kmer_list;
    vector<Kmer_Occurence_Pair> kmer_list;
    
    Kmer_counter_map_iterator it;
    
    cerr << "Kcounter hash size: " << num_kmers << endl;
    
    unsigned long count = 0;

    for (it = _kmer_counter.begin(); it != _kmer_counter.end(); it++) {

        if (it->second > 0) {

            Kmer_Occurence_Pair kp (it->first, it->second);
        
            kmer_list.push_back(kp);
            count++;
        }
    }
    
    cerr << "Processed " << count << " non-zero abundance kmers in kcounter.\n";
    
    if (IRKE_COMMON::__DEVEL_no_kmer_sort) {
        
        cerr << "-Not sorting list of kmers, given parallel mode in effect." << endl;
    
        return(kmer_list);
    }
    
    cerr << "Sorting " << num_kmers << " kmers ...";

    Sort_kmer_by_count_desc sorter (this);
   
    sort(kmer_list.begin(), kmer_list.end(), sorter);
    
    end = time (NULL);
    
    long time_spent = end - start;
    
    cerr << "Done sorting " << num_kmers << " kmers, taking " << time_spent << " seconds." << endl;
    
    return(kmer_list);
}



/******************/
/*  Kmer_visitor */
/****************/

Kmer_visitor::Kmer_visitor(unsigned int kmer_length, bool is_ds)
{
	_kmer_length = kmer_length;
	_DS_MODE = is_ds;
}

void Kmer_visitor::add (kmer_int_type_t kmer)
{
	if (_DS_MODE) {
		kmer = get_DS_kmer_val(kmer, _kmer_length);
	}
	_set.insert(kmer);
}

bool Kmer_visitor::exists (kmer_int_type_t kmer)
{
	if (_DS_MODE) {
        kmer = get_DS_kmer_val(kmer, _kmer_length);
    }
	return (_set.find(kmer) != _set.end());
}

void Kmer_visitor::erase (kmer_int_type_t kmer)
{
	if (_DS_MODE)
	{
		kmer = get_DS_kmer_val(kmer, _kmer_length);
	}
	_set.erase(kmer);
}

void Kmer_visitor::clear()
{
	_set.clear();
}

unsigned int Kmer_visitor::size()
{
	return(_set.size());
}



/******************/
/*  Sort Classes */
/****************/

Sort_kmer_by_count_desc::Sort_kmer_by_count_desc(KmerCounter *kcounter)
{
	this->kcounter = kcounter;
}

bool Sort_kmer_by_count_desc::operator() (const Kmer_counter_map_iterator& i, const Kmer_counter_map_iterator& j) {

#ifdef _DEBUG
    if (i->second == j->second)
        return (i->first > j->first);
	else
#endif
        return (i->second > j->second);

}

bool Sort_kmer_by_count_desc::operator() (const Kmer_Occurence_Pair& i, const Kmer_Occurence_Pair& j) {

#ifdef _DEBUG
    if (i.second == j.second)
        return (i.first > j.first);
	else
#endif
        return (i.second > j.second);

}

bool Sort_kmer_by_count_desc::operator() (const kmer_int_type_t& val_i, const kmer_int_type_t& val_j) {

    unsigned int count_i = kcounter->get_kmer_count(val_i);
    unsigned int count_j = kcounter->get_kmer_count(val_j);
    
#ifdef _DEBUG
    if (count_i == count_j)
        return (val_i > val_j);
	else
#endif
        return (count_i > count_j);
    
}

bool Sort_kmer_by_count_desc::operator() (const string& i, const string& j) {
    
	kmer_int_type_t val_i = kcounter->get_kmer_intval(i);
	kmer_int_type_t val_j = kcounter->get_kmer_intval(j);
	unsigned int count_i = kcounter->get_kmer_count(val_i);
    unsigned int count_j = kcounter->get_kmer_count(val_j);

#ifdef _DEBUG
    if (count_i == count_j)
        return (val_i > val_j);
	else
#endif
        return (count_i > count_j);
}
