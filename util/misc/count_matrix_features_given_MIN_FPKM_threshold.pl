#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 fpkm.matrix\n\n";


my $matrix_file = $ARGV[0] or die $usage;

open (my $fh, $matrix_file) or die $!;
my $header = <$fh>;

my @fpkms;
while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    shift @x; # gene accession
    my $max_fpkm = shift @x;
    while (@x) {
        my $fpkm = shift @x;
        if ($fpkm > $max_fpkm) {
            $max_fpkm = $fpkm;
        }
    }
    push (@fpkms, $max_fpkm);
}

@fpkms = reverse sort {$a<=>$b} @fpkms;

my $min_fpkm_thresh = int($fpkms[0] + 0.5);
my $num_features = 1;

print "neg_min_fpkm\tnum_features\n";

shift @fpkms;
while (@fpkms) {

    my $fpkm = shift @fpkms;
    $fpkm = int($fpkm+0.5);

    if ($fpkm < $min_fpkm_thresh) {
        print "" . (-1*$min_fpkm_thresh) . "\t$num_features\n";
        $min_fpkm_thresh = $fpkm;

    }
    $num_features++;
}

print "$min_fpkm_thresh\t$num_features\n";

exit(0);
