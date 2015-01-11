#!/usr/bin/env perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Fasta_reader;


my $usage = "usage: $0 seqs.fasta\n\n";

my $file = $ARGV[0] or die $usage;

my %seq_to_header;

my $fasta_reader = new Fasta_reader($file);
while (my $seq_obj = $fasta_reader->next()) {

    my $sequence = $seq_obj->get_sequence();
    my $header = $seq_obj->get_header();
    
    if (exists $seq_to_header{$sequence}) {
        $seq_to_header{$sequence} .= "\t$header";
    }
    else {
        $seq_to_header{$sequence} = $header;
    }
}

foreach my $sequence (keys %seq_to_header) {
    my $header = $seq_to_header{$sequence};

    print ">$header\n$sequence\n";
}

exit(0);


