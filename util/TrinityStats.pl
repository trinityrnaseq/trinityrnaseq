#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../PerlLib");
use Fasta_reader;
use BHStats;

my $usage = "\n\nusage: $0 transcripts.fasta\n\n";

my $fasta_file = $ARGV[0] or die $usage;

main: {

    my $fasta_reader = new Fasta_reader($fasta_file);
    
    my @all_seq_lengths;

    my $number_transcripts = 0;
    
    my $num_GC = 0;
        
    my %component_to_longest_isoform;
    
    my $tot_seq_len = 0;

    my $missing_gene_ids_flag = 0;
    
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        
        $number_transcripts++;
         my $comp_id = $acc;


        if (! $missing_gene_ids_flag) {
            
            if ($acc =~ /^(.*c\d+_g\d+)/) {
                $comp_id = $1;
            }
            elsif ($acc =~ /^(.*comp\d+_c\d+)/) {
                $comp_id = $1;
            }
            else {
                print STDERR "Error, cannot decipher gene identifier from acc: $acc";
                $missing_gene_ids_flag = 1;
            }
        }
        
        my $sequence = $seq_obj->get_sequence();

        my $seq_len = length($sequence);

        $tot_seq_len += $seq_len;

        if ( (! exists $component_to_longest_isoform{$comp_id})
             || 
             $component_to_longest_isoform{$comp_id} < $seq_len) {
            
            $component_to_longest_isoform{$comp_id} = $seq_len;
        }
        
        push (@all_seq_lengths, $seq_len);
                
        while ($sequence =~ /[gc]/ig) {
            $num_GC++;
        }

    }

    print "\n\n";
    print "################################\n";
    print "## Counts of transcripts, etc.\n";
    print "################################\n";

    print "Total trinity 'genes':\t" . scalar(keys %component_to_longest_isoform) . "\n";        
    print "Total trinity transcripts:\t" . $number_transcripts . "\n";

    
    my $pct_gc = sprintf("%.2f", $num_GC / $tot_seq_len * 100);
    

    print "Percent GC: $pct_gc\n\n"; 
    
    print "########################################\n";
    print "Stats based on ALL transcript contigs:\n";
    print "########################################\n\n";

    &report_stats(@all_seq_lengths);
    print "\n\n";
    
    if ($missing_gene_ids_flag) {
        print " - note: not reporting gene-based longest isoform info since couldn't parse Trinity accession info.\n";
    }
    else {
        print "#####################################################\n";
        print "## Stats based on ONLY LONGEST ISOFORM per 'GENE':\n";
        print "#####################################################\n\n";
        
        &report_stats(values %component_to_longest_isoform);
            
        print "\n\n\n";
    }
    
    
    exit(0);
}

####
sub report_stats {
    my (@seq_lengths) = @_;

    @seq_lengths = reverse sort {$a<=>$b} @seq_lengths;
    
    my $cum_seq_len = 0;
    foreach my $len (@seq_lengths) {
        $cum_seq_len += $len;
    }
    
    for (my $i = 10; $i <= 50; $i += 10) {
        my $cum_len_needed = $cum_seq_len * $i/100;
        my $N_val = &get_contigNvalue($cum_len_needed, \@seq_lengths);
        print "\tContig N$i: $N_val\n";
    }
    print "\n";
    

    my $median_len = &BHStats::median(@seq_lengths);
    print "\tMedian contig length: $median_len\n";

    my $avg_len = sprintf("%.2f", &BHStats::avg(@seq_lengths));
    print "\tAverage contig: $avg_len\n";

    
        
    print "\tTotal assembled bases: $cum_seq_len\n";

    return;
}



sub get_contigNvalue {
    my ($cum_len_needed, $seq_lengths_aref) = @_;
    
    my $partial_sum_len = 0;
    foreach my $len (@$seq_lengths_aref) {
        $partial_sum_len += $len;

        if ($partial_sum_len >= $cum_len_needed) {
            return($len);
        }
    }


    return -1; # shouldn't happen.
}
