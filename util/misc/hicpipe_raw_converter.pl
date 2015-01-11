#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0  bowtie.raw  \n\n";

my $bowtie_raw = $ARGV[0] or die $usage;

open (my $fh, $bowtie_raw) or die "Error, cannot open file $bowtie_raw";

my $header = <$fh>;
print $header;

while (<$fh>) {
    if (/random/) { next; }
    if (/chrUn/) { next; }
    s/chr//g;
    print;
}

exit(0);

