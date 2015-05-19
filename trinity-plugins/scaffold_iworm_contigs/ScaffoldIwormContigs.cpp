#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <stdexcept>

#ifndef WIN32
#define HTS
#endif

#ifdef HTS
#include "htslib/sam.h"
#endif

using namespace std;

const char *my_endl = "\n";

struct map_key
{
	template <typename T>
	typename T::first_type operator()(T keyValuePair) const
	{
		return keyValuePair.first;
	}
};

void examine_frags(const set<string>& iworm_left, const set<string>& iworm_right, map<pair<string, string>, int>& paired_iworm_contigs);

ofstream SAM_OFH("scaffolding_entries.sam");

template <typename T>
struct sort_pred {
	map<T, int> &_hash;
	sort_pred(map<T, int>& h) : _hash(h) {}
	bool operator()(const T &left, const T &right) {
		return _hash[left] > _hash[right];
	}
};

template <typename T>
void copy_to_out(map<T, int>& paired_iworm_contigs, map<string, int>& iworm_acc_to_fasta_index)
{
	vector<pair<string,string> > outvec;
	transform(paired_iworm_contigs.begin(), paired_iworm_contigs.end(), back_inserter(outvec), map_key());
	sort(outvec.begin(), outvec.end(), sort_pred<pair<string, string> >(paired_iworm_contigs));
	for (vector<pair<string, string> >::const_iterator a = outvec.begin(); a != outvec.end(); ++a)
	{
		int index_A = iworm_acc_to_fasta_index[a->first];
		int index_B = iworm_acc_to_fasta_index[a->second];

		cout << a->first << '\t' << index_A << '\t' << a->second << '\t' << index_B << '\t' << paired_iworm_contigs[*a] << my_endl;
	}
}

void index_iworm(string iworm_fasta, map<string, int>& aMap)
{
	int counter = 0;

	ifstream fh(iworm_fasta.c_str());
	while (fh)
	{
		string line;
		fh >> line;
		if (line[0] == '>')
			aMap[line.substr(1, string::npos)] = counter++;
	}
}

void process_line(string read_acc, string iworm_acc, set<string>& iworm_left, 
	set<string>& iworm_right, string& prev_core_acc,
	map<pair<string, string>, int>& paired_iworm_contigs,
	int& num_unrecognized_iworm_contig_names, int& num_warnings)
{
	string parsable(iworm_acc);
	if (parsable[0] == 'a') parsable[0] = ' ';
	replace(parsable.begin(), parsable.end(), ';', ' ');
	istringstream iss2(parsable);
	int a1, a2;
	iss2 >> a1 >> a2;
	if (iss2.fail())
	{
		num_unrecognized_iworm_contig_names++;
		if (num_unrecognized_iworm_contig_names <= 10) {
			cerr << "warning, inchworm contig (" << iworm_acc << ") isn't recognized. Ignoring entry: [[" << read_acc << "]]";
		}
		if (num_unrecognized_iworm_contig_names == 11) {
			cerr << "warning, too many unrecognized inchworm contig names.  Will report summary of counts later.";
		}
		return;
	}

	int s = read_acc.rfind('/');
	string core_acc = read_acc.substr(0, s);
	char frag_end = read_acc[s + 1] - '1';

	if (s == string::npos || (frag_end != 0 && frag_end != 1))
	{
		// must have mixed in a single read with the pairs...
		if (num_warnings <= 10) {
			cerr << "warning, ignoring read: " << read_acc << "since cannot decipher if /1 or /2 of a pair.";
		}
		else if (num_warnings == 11) {
			cerr << "number of read warnings exceeded 10.  Turning off warning messages from here out.\n";
		}

		num_warnings++;
		return;
	}

	if (core_acc != prev_core_acc)
	{
		examine_frags(iworm_left, iworm_right, paired_iworm_contigs);
		iworm_left.clear();
		iworm_right.clear();
	}

	if (frag_end == 0)
		iworm_left.insert(iworm_acc);
	if (frag_end == 1)
		iworm_right.insert(iworm_acc);

	prev_core_acc = core_acc;
}

int main(int argc, char* argv[])
{
	string usage = "\n\nusage: ScaffoldIwormContigs nameSorted.sam iworm_fasta_file\n\n";
	if (argc != 3)
	{
		cerr << usage;
		getchar();
		return -1;
	}

	string name_sorted_sam_file(argv[1]);
	map<string, int> iworm_acc_to_fasta_index;
	index_iworm(argv[2], iworm_acc_to_fasta_index);

	int num_unrecognized_iworm_contig_names = 0, num_warnings = 0;
	string prev_core_acc;
	map<pair<string,string>, int> paired_iworm_contigs;
	set<string> iworm_left;
	set<string> iworm_right;

#ifdef HTS
	samFile *in = sam_open(name_sorted_sam_file.c_str(), "r");
	bam_hdr_t *h = sam_hdr_read(in);
	bam1_t *b = bam_init1();

	while (sam_read1(in, h, b) > 0)
	{
		string read_acc(bam_get_qname(b));
		if (b->core.tid >= 0)
		{
			string iworm_acc(h->target_name[b->core.tid]);
			process_line(read_acc, iworm_acc, iworm_left, iworm_right, prev_core_acc, paired_iworm_contigs, num_unrecognized_iworm_contig_names, num_warnings);
		}
	}

	int r = sam_close(in);
	if (r < 0) {
		throw runtime_error("Error closing input");
	}
#else
	ifstream fh2(name_sorted_sam_file.c_str());
	string line;
	while (std::getline(fh2, line))
	{
		istringstream iss(line);
		string read_acc, tmp, iworm_acc;
		iss >> read_acc >> tmp >> iworm_acc;
		if (iworm_acc.empty())
			continue;
		if (iworm_acc == "*")
			continue;

		process_line(read_acc, iworm_acc, iworm_left, iworm_right, prev_core_acc, paired_iworm_contigs, num_unrecognized_iworm_contig_names, num_warnings);
	}
#endif
	// get last one
	examine_frags(iworm_left, iworm_right, paired_iworm_contigs);


	if (num_warnings) {
		cerr << "WARNING: note there were " << num_warnings << " reads that could not be deciphered as being /1 or /2 of a PE fragment.  Hopefully, these were SE reads that should have been ignored. Otherwise, please research this further.\n\n";
	}
	if (num_unrecognized_iworm_contig_names) {
		cerr <<  "WARNING: note, there were " << num_unrecognized_iworm_contig_names << " inchworm contig names in the SAM file that were ignored due to the inchworm contig accession not being recognized.\n";
	}

	copy_to_out(paired_iworm_contigs, iworm_acc_to_fasta_index);

#ifdef _DEBUG
	getchar();
#endif
	return 0;
}

void examine_frags(const set<string>& iworm_left, const set<string>& iworm_right, map<pair<string, string>, int>& paired_iworm_contigs)
{
	if (iworm_left.size() == 0 || iworm_right.size() == 0)
		return;	// no pairs

	if (iworm_left.size() == 1 && iworm_right.size() == 1)
	{
		string k1 = *iworm_left.begin();
		string k2 = *iworm_right.begin();
		if (k1 != k2)
		{
			if (k2 < k1) swap(k1, k2);

			pair<string, string> key(k1, k2);
			if (paired_iworm_contigs.find(key) == paired_iworm_contigs.end())
				paired_iworm_contigs[key] = 1;
			else
				paired_iworm_contigs[key]++;

			SAM_OFH << k1 << k2 << endl;
		}
	}
}

