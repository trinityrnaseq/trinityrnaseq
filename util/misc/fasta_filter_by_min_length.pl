#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;


my $usage = "usage: $0 fastaFile min_length\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $min_length = $ARGV[1] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);

while (my $seqobj = $fasta_reader->next()) {
    my $fasta_entry = $seqobj->get_FASTA_format();
    my $sequence = $seqobj->get_sequence();
    unless (length($sequence) >= $min_length) {
        next;
    }

    print $fasta_entry;
}


exit(0);

