#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 trans_gff3.list\n\n";

my $trans_gff3_files = $ARGV[0] or die $usage;

main: {

    open (my $fh, $trans_gff3_files) or die $!;

    while (<$fh>) {
        chomp;
        my $filename = $_;

        my $cmd = "transcript_gff3_to_bed.pl $filename > $filename.bed";
        
        print "$cmd\n";
    }

    close $fh;

    exit(0);
}


