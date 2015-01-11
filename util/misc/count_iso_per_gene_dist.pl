#!/usr/bin/env perl

use strict;
use warnings;


my %counter;

while (<>) {
    if (/^>/) {
        if (/^>(.*c\d+_g\d+)/) {
            $counter{$1}++;
        }
        elsif (/^>(.*comp\d+_c\d+)/) {
            # older format
            $counter{$1}++;
        }
        else {
            die "Error, dont recognize formatting of accession in line: $_";
        }
    }
}

my %count_counter;
for my $val (values %counter) {
    $count_counter{$val}++;
}

print join("\t", "#iso_per_gene", "num_genes") . "\n";

for my $count (sort {$a<=>$b} keys %count_counter) {
    my $val = $count_counter{$count};

    print join("\t", $count, $val) . "\n";
}

exit(0);
    
