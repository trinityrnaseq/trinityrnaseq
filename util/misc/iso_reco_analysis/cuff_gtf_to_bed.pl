#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 cuff_gtf.list\n\n";

my $trans_gtf_files = $ARGV[0] or die $usage;

main: {

    open (my $fh, $trans_gtf_files) or die $!;

    while (<$fh>) {
        chomp;
        my $filename = $_;

        my $cmd = "cufflinks_gtf_to_bed.pl $filename > $filename.bed";
        
        print "$cmd\n";
    }
    
    close $fh;

    exit(0);
}


