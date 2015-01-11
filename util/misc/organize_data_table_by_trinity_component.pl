#!/usr/bin/env perl

use strict;
use warnings;


my %data;
while (<>) {
    if (/^\#/) {
        print;
        next;
    }
    unless (/\w/) { next; }

    my $line = $_;
    if (/^(comp\d+_c\d+)/) {
        my $comp = $1;
        $data{$comp} .= $line;
    }
    else {
        die "Error, cannot decode component identity from $line";
    }
}

foreach my $component (keys %data) {
    
    print $data{$component} . "\n";

}

exit(0);

