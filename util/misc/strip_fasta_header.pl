#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 file.fasta\n\n";

my $fasta_file = $ARGV[0] or die $usage;

main: {


    my $fasta_reader = new Fasta_reader($fasta_file);
    
    while (my $seq_obj = $fasta_reader->next()) {

        my $sequence = $seq_obj->get_sequence();

        print ">\n$sequence\n";
    }


    exit(0);
}
