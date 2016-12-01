#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../../PerlLib");
use PSL_parser;
use DelimParser;

my $usage = "\n\n\tusage: $0 blat_output.pslx\n\n";

my $blat_output = $ARGV[0] or die $usage;

main: {

    my @column_fields = ("target", "query", "per_id", "per_gap", "score",
                         "target_lend", "target_rend", "query_lend", "query_rend",
                         "strand", "per_len");
    
    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", \@column_fields);
    
    
    my $psl_parser = new PSL_parser($blat_output);
    
	while (my $psl_entry = $psl_parser->get_next()) {
		
        #print $psl_entry->toString();
        
		my $trans_assembly = $psl_entry->get_Q_name();
		my $gene_id = $psl_entry->get_T_name();
		my $per_id = $psl_entry->get_per_id();
		
		my $num_gap_opens = $psl_entry->get_T_gap_count() + $psl_entry->get_Q_gap_count();
		my $num_gap_bases = $psl_entry->get_T_gap_bases() + $psl_entry->get_Q_gap_bases();
		my $num_matches = $psl_entry->get_match_count();
		my $num_mismatches = $psl_entry->get_mismatch_count();
		
        
		my ($trans_end5, $trans_end3) = $psl_entry->get_Q_span();
		my ($gene_end5, $gene_end3) = $psl_entry->get_T_span();
		
		my $orient = $psl_entry->get_strand();
		
		
		my $gene_seq_len = $psl_entry->get_T_size();
		
		my $percent_gapped = ($num_gap_bases) / ($num_matches + $num_mismatches) * 100;
		
		my ($gene_lend, $gene_rend) = sort {$a<=>$b} ($gene_end5, $gene_end3);
        
		my ($trans_lend, $trans_rend) = sort {$a<=>$b} ($trans_end5, $trans_end3);
        
		my $delta = ($gene_lend - 1) + ($gene_seq_len - $gene_rend);
		
        my $pct_len = 100 - ($delta/$gene_seq_len * 100);
        $pct_len = sprintf("%.2f", $pct_len);
        
        ## score the alignment:
        my $score = (5 * $num_matches) - (4 * $num_mismatches) - (20 * $num_gap_opens) - log($num_gap_bases + 1);
        
        
        my %row = ( 
            
            target => $gene_id,
            query => $trans_assembly,
            
            per_id => $per_id,
            
            score => sprintf("%.2f", $score),
            
            target_lend => $gene_lend,
            target_rend => $gene_rend,
            
            query_lend => $trans_lend,
            query_rend => $trans_rend,
            
            strand => $orient,

            per_gap => sprintf("%.2f", $percent_gapped),
            
            per_len => $pct_len,
            
            #psl => $psl_entry,
            
            );
        
        $tab_writer->write_row(\%row);
    }
     
    
}

