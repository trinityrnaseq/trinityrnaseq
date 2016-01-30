#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;


my $usage = "\n\nusage: $0 transcripts.fasta\n\n";

my $fasta_file = $ARGV[0] or die $usage;

main: {

    my $fasta_reader = new Fasta_reader($fasta_file);
    
    my @seq_lengths;
    my $cum_seq_len = 0;
    
    while (my $seq_obj = $fasta_reader->next()) {

        my $sequence = $seq_obj->get_sequence();

        my $seq_len = length($sequence);

        $cum_seq_len += $seq_len;
        push (@seq_lengths, $seq_len);
    }

    @seq_lengths = reverse sort {$a<=>$b} @seq_lengths;

    my $half_cum_len = $cum_seq_len / 2;

    my $partial_sum_len = 0;
    foreach my $len (@seq_lengths) {
        $partial_sum_len += $len;

        if ($partial_sum_len >= $half_cum_len) {
            print "N50: $len\n";
            last;
        }
    }

    exit(0);
}

