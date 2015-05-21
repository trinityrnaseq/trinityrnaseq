#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include "error_checker.h"

#ifndef WIN32
#define HTS
#endif

#ifdef HTS
#include "htslib/sam.h"
#endif

#ifdef _DEBUG
ofstream SAM_OFH("scaffolding_entries.sam");
#endif

using namespace std;


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

#ifdef _DEBUG
			SAM_OFH << k1 << k2 << endl;
#endif

		}
	}
}

struct map_key
{
	template <typename T>
	typename T::first_type operator()(T keyValuePair) const
	{
		return keyValuePair.first;
	}
};

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
	vector<T> outvec;
	const char *my_endl = "\n";

	outvec.reserve(paired_iworm_contigs.size());
	transform(paired_iworm_contigs.begin(), paired_iworm_contigs.end(), back_inserter(outvec), map_key());
	sort(outvec.begin(), outvec.end(), sort_pred<T >(paired_iworm_contigs));
	for (typename vector<T>::const_iterator a = outvec.begin(); a != outvec.end(); ++a)
	{
		const int index_A = iworm_acc_to_fasta_index[a->first];
		const int index_B = iworm_acc_to_fasta_index[a->second];

		cout << a->first << '\t' << index_A << '\t';
		cout << a->second << '\t' << index_B << '\t';
		cout << paired_iworm_contigs[*a] << my_endl;
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

void validate_iworm(string iworm_acc)
{
	pair<int, int> a;

	if (iworm_acc[0] == 'a') iworm_acc[0] = ' ';
	replace(iworm_acc.begin(), iworm_acc.end(), ';', ' ');
	istringstream iss2(iworm_acc);
	iss2 >> a.first >> a.second;
	if (iss2.fail())
		throw bad_name_exception();
}

pair<string, int> parse_read_name(string read_acc)
{
	pair<string, int> result;
	size_t s = read_acc.rfind('/');
	result.first = read_acc.substr(0, s);
	result.second = read_acc[s + 1] - '1';

	if (s == string::npos || (result.second != 0 && result.second != 1))
		throw single_read_exception();

	return result;
}

void process_line(string read_acc, string iworm_acc, set<string>& iworm_left, 
	set<string>& iworm_right, string& prev_core_acc,
	map<pair<string, string>, int>& paired_iworm_contigs,
	error_checker& checker)
{
	try
	{
		validate_iworm(iworm_acc);

		pair<string, int> read = parse_read_name(read_acc);
		const string core_acc = read.first;
		int frag_end = read.second;

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
	catch (const bad_name_exception&)
	{
		checker.add_bad_name(iworm_acc, read_acc);
	}
	catch (const single_read_exception&)
	{
		checker.add_single_read(read_acc);
	}
}

int main(int argc, char* argv[])
{
	string usage = "\n\nusage: ScaffoldIwormContigs nameSorted.sam iworm_fasta_file\n\n";
	if (argc != 3)
	{
		cerr << usage;
		return -1;
	}

	string name_sorted_sam_file(argv[1]);
	map<string, int> iworm_acc_to_fasta_index;
	index_iworm(argv[2], iworm_acc_to_fasta_index);

	string prev_core_acc;
	map<pair<string,string>, int> paired_iworm_contigs;
	set<string> iworm_left;
	set<string> iworm_right;
	error_checker errors;

	samFile *in = sam_open(name_sorted_sam_file.c_str(), "r");
	bam_hdr_t *h = sam_hdr_read(in);
	bam1_t *b = bam_init1();

	while (sam_read1(in, h, b) > 0)
	{
		string read_acc(bam_get_qname(b));
		if (b->core.tid >= 0)
		{
			string iworm_acc(h->target_name[b->core.tid]);
			process_line(read_acc, iworm_acc, iworm_left, iworm_right, prev_core_acc, paired_iworm_contigs, errors);
		}
	}

	int r = sam_close(in);
	if (r < 0) {
		throw runtime_error("Error closing input");
	}

	// get last one
	examine_frags(iworm_left, iworm_right, paired_iworm_contigs);

	errors.report();

	copy_to_out(paired_iworm_contigs, iworm_acc_to_fasta_index);

	return 0;
}


