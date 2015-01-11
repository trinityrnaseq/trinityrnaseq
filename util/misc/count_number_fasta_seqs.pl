#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
    my $filename = $_;
    chomp $filename;

    my $count = 0;
    
    open (my $fh, $filename) or die "Error, cannot open file $filename";
    while (<$fh>) {
        if (/^>/) {
            $count++;
        }
    }
    print "$count\t$filename\n";
    close $fh;
    
}


exit(0);

