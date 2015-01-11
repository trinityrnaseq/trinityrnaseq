#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
    if (/Welding: >(a\d+;\d+).*to >(a\d+;\d+)/) {
        my $from = $1;
        my $to = $2;
        print join("\t", sort($from, $to)) . "\n";
    }
    elsif (/SCAFFOLD_ACCEPT:\s+(a\d+;\d+)\s+\d+\s+(a\d+;\d+)/) {
        my $from = $1;
        my $to = $2;
        print join("\t", sort($from, $to)) . "\n";
    }
}

exit(0);

