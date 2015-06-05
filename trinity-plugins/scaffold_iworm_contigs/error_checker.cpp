#include <string>
#include <iostream>
#include "error_checker.h"

using namespace std;

void error_checker::add_bad_name(string contig, string entry)
{
	num_unrecognized_iworm_contig_names++;
	if (num_unrecognized_iworm_contig_names <= 10) {
		cerr << "warning, inchworm contig (" << contig << ") isn't recognized. Ignoring entry: [[" << entry << "]]";
	}
	if (num_unrecognized_iworm_contig_names == 11) {
		cerr << "warning, too many unrecognized inchworm contig names.  Will report summary of counts later.";
	}
}

void error_checker::add_single_read(string read)
{
	// must have mixed in a single read with the pairs...
	if (num_warnings <= 10) {
		cerr << "warning, ignoring read: " << read << "since cannot decipher if /1 or /2 of a pair.";
	}
	else if (num_warnings == 11) {
		cerr << "number of read warnings exceeded 10.  Turning off warning messages from here out.\n";
	}

	num_warnings++;
}

void error_checker::report() const
{
	if (num_warnings) {
		cerr << "WARNING: note there were " << num_warnings << " reads that could not be deciphered as being /1 or /2 of a PE fragment.  Hopefully, these were SE reads that should have been ignored. Otherwise, please research this further.\n\n";
	}
	if (num_unrecognized_iworm_contig_names) {
		cerr << "WARNING: note, there were " << num_unrecognized_iworm_contig_names << " inchworm contig names in the SAM file that were ignored due to the inchworm contig accession not being recognized.\n";
	}
}

