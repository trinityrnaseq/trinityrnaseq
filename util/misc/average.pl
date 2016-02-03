#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use BHStats;

my $count = 0;
my $sum = 0;
my @values;

while (<STDIN>) {
    chomp;
    $sum += $_;
    $count++;
    push (@values, $_);
}

@values = sort {$a<=>$b} @values;

if ($count) {
    my $average = ($sum/$count);
    my $median_pos = int ($count/2);
    print "\n"; 
    print "MIN: " . $values[0] . "\n";
    print "MAX: " . $values[$#values] . "\n";
    print "Sum: $sum\n";
    printf ("Average: %.2f\n", $average);
    print "Median: " . BHStats::median(@values) . "\n";
    my $stdev = BHStats::stDev(@values);
    printf ("stDev from Average: %.2f\n", $stdev);

    my $geoMean = &BHStats::geometric_mean(@values);
    if ($geoMean) {
	printf ("geoMean: %.2f\n", $geoMean);
    }
}

