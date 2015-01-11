#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 subcluster_A subcluster_B ... >  merged.subclusters\n\n";

my @files = @ARGV;
unless (@files) {
    die $usage;
}

my $printed_header_flag = 0;

foreach my $file (@files) {
    
    open (my $fh, $file) or die "Error, cannot open file $file";
    my $header = <$fh>;
    
    print $header unless $printed_header_flag;
    $printed_header_flag = 1;

    while (<$fh>) {
        print;
    }
    close $fh;
}

exit(0);

