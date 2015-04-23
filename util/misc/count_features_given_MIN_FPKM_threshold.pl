#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 expr.RSEM\n\n";

my $expr_file = $ARGV[0] or die $usage;

open (my $fh, $expr_file) or die $!;
my $header = <$fh>;

my @fpkms;
while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $fpkm = $x[6];
    push (@fpkms, $fpkm);
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

print "" . (-1*$min_fpkm_thresh) . "\t$num_features\n";

exit(0);
