#!/usr/bin/env perl

use strict;
use warnings;

use Fasta_retriever;
use List::Util qw(shuffle);

my $usage = "usage: $0 file.fasta\n\n";

my $file = $ARGV[0] or die $usage;

main: {

    my @accs;
    open (my $fh, $file) or die $!;
    while (<$fh>) {
        if (/^>(\S+)/) {
            my $acc = $1;
            push (@accs, $acc);
        }
    }
    close $fh;


    @accs = shuffle @accs;

    my $fasta_retriever = new Fasta_retriever($file);

    foreach my $acc (@accs) {
        
        my $seq = $fasta_retriever->get_seq($acc);

        print "$acc\t$seq\n";
    }

    exit(0);
}

