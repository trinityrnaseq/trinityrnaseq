#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use Process_cmd;

my $usage = "\n\n\tusage: $0 transcript_aligned.bam [out_prefix='ss_analysis']\n\n";

my $bam_file = $ARGV[0] or die $usage;
my $out_prefix = $ARGV[1] || "ss_analysis";

main: {

    my %transcript_to_orients;

    my $sam_reader = new SAM_reader($bam_file);
    print STDERR "-parsing file: $bam_file\n";
    while (my $sam_entry = $sam_reader->get_next()) {

        if ($sam_entry->is_proper_pair() && $sam_entry->is_first_in_pair()) {

            my $orient = $sam_entry->get_query_strand();
            my $trans_name = $sam_entry->get_scaffold_name();

            #print STDERR "$trans_name\t$orient\n";
            
            $transcript_to_orients{$trans_name}->{$orient}++;
        }
    }
    print STDERR "-done parsing file, examining orientations of reads.\n";

    ## sum them up.
    my @transcripts = keys %transcript_to_orients;
    foreach my $transcript (@transcripts) {
        
        my $orient_plus = $transcript_to_orients{$transcript}->{'+'} || 0;
        my $orient_minus = $transcript_to_orients{$transcript}->{'-'} || 0;
        

        my $total_reads = $orient_plus + $orient_minus;
        
        $transcript_to_orients{$transcript}->{'transcript'} = $transcript;
        $transcript_to_orients{$transcript}->{'total_reads'} = $total_reads;
    }

    ####
    my @structs = values %transcript_to_orients;
    @structs = reverse sort {$a->{total_reads}<=>$b->{total_reads}} @structs; 
    
    # header
    open (my $ofh, ">$out_prefix.dat") or die "Error, cannot write to $out_prefix.dat";
    
    print $ofh join("\t", "#transcript", "plus_strand_1stReads", "minus_strand_1stReads", "total_reads", "diff_ratio") . "\n";
    
    foreach my $struct (@structs) {
        
        my $struct = shift @structs;
        

        my ($plus, $minus, $total_reads, $transcript) = ($struct->{'+'},
                                                         $struct->{'-'},
                                                         $struct->{'total_reads'},
                                                         $struct->{'transcript'});
        unless ($plus) { $plus = 0; }
        unless ($minus) { $minus = 0; }

        my $diff_proportion = sprintf("%.3f", ($plus - $minus) / $total_reads);
        
        print $ofh join("\t", $transcript, $plus, $minus, $total_reads, $diff_proportion) . "\n";
        
    }
    close $ofh;


    ## plot it.
    my $cmd  = "$FindBin::Bin/plot_strand_specificity_dist_by_quantile.Rscript $out_prefix.dat";
    &process_cmd($cmd);
    
    
    exit(0);
}


          
