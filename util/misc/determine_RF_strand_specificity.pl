#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\n\n\tusage: $0 transcript_aligned.bam top_percent=10\n\n";

my $bam_file = $ARGV[0] or die $usage;
my $top_percent = $ARGV[1] || 10;

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
    
    my $num_transcripts = scalar(@structs);
    my $max_index = int($top_percent/100 * $num_transcripts);

    # header
    print join("\t", "#transcript", "minus_strand_1stReads", "plus_strand_1stReads", "total_reads", "pct_RF_SS") . "\n";
    
    for (0..$max_index) {
        my $struct = shift @structs;
        

        my ($plus, $minus, $total_reads, $transcript) = ($struct->{'+'},
                                                         $struct->{'-'},
                                                         $struct->{'total_reads'},
                                                         $struct->{'transcript'});

        # assuming RF
        my $pct_total = $minus / $total_reads * 100;
        
        print join("\t", $transcript, $minus, $plus, $total_reads, sprintf("%.1f", $pct_total)) . "\n";

    }
    
    exit(0);
}


          
