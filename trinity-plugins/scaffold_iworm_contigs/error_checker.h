#ifndef ERROR_CHECKER_H_B91F75504340473B989F054E04E58BE7
#define ERROR_CHECKER_H_B91F75504340473B989F054E04E58BE7

#include <string>

class bad_name_exception : std::exception
{
};

class single_read_exception : std::exception
{
};

struct error_checker
{
	void add_bad_name(std::string contig, std::string entry);
	void add_single_read(std::string read);
	void report() const;
private:
	int num_unrecognized_iworm_contig_names;
	int num_warnings;
};

#endif
