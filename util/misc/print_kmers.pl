#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 file.fa [kmer_length=25]\n\n";

my $fa_file = $ARGV[0] or die $usage;
my $kmer_length = $ARGV[1] || 25;

main: {

    my $fasta_reader = new Fasta_reader($fa_file);
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $sequence = $seq_obj->get_sequence();
        
        for (my $i = 0; $i < length($sequence) - $kmer_length + 1; $i++) {
            
            my $kmer = substr($sequence, $i, $kmer_length);

            print "$kmer\n";
        }
    }

    exit(0);
}
