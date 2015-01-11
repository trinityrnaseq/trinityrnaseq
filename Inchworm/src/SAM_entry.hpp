#ifndef __SAM_ENTRY__

#define __SAM_ENTRY__

#include <string>
#include <vector>
#include <utility>

typedef struct {
    unsigned long genome_lend;
    unsigned long genome_rend;
    
    unsigned long query_lend;
    unsigned long query_rend;
    
} alignment_segment;



using namespace std;

class SAM_entry {
    
public:
    
    SAM_entry(string sam_line);
    
    string get_sam_string();
    
    string toString();
    
    string get_read_name();
    
    string get_full_read_name();
    
    string get_scaffold_name();
    
    unsigned long get_scaffold_position();
    
    void set_scaffold_position(unsigned long position);
    
    unsigned int get_mapping_quality();
    
    string get_cigar_alignment();
    
    void set_cigar_alignment(string cigar_string);
    
    vector<alignment_segment> get_alignment_coords();
    
    string get_mate_scaffold_name();
    
    unsigned long get_mate_scaffold_position();
    
    string get_sequence();
    
    string get_quality_scores();
    
    unsigned int get_flag();
    
    bool is_paired();
    
    bool is_proper_pair();
    
    bool is_query_unmapped();
    
    bool is_mate_unmapped();
    
    char get_query_strand();
    
    char get_query_transcribed_strand(); // takes into account strand-specific sequencing and left/right pair info assuming our dUTP protocol.
    
    char get_mate_strand();
    
    bool is_first_in_pair();
    
    bool is_second_in_pair();
    
private:
  string line;
  vector<string> tokens;
  

};


#endif

