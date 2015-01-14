#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 gene_gff3.list\n\n";

my $gene_gff3_files = $ARGV[0] or die $usage;

main: {

    open (my $fh, $gene_gff3_files) or die $!;

    while (<$fh>) {
        chomp;
        my $filename = $_;

        my $cmd = "gene_gff3_to_bed.pl $filename > $filename.bed";

        print "$cmd\n";
    }

    close $fh;

    exit(0);
}


