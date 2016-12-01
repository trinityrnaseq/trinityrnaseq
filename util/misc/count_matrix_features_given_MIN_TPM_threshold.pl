#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 tpm.matrix\n\n";


my $matrix_file = $ARGV[0] or die $usage;

open (my $fh, $matrix_file) or die $!;
my $header = <$fh>;

my @tpms;
while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    shift @x; # gene accession
    my $max_tpm = shift @x;
    while (@x) {
        my $tpm = shift @x;
        if ($tpm > $max_tpm) {
            $max_tpm = $tpm;
        }
    }
    push (@tpms, $max_tpm);
}

@tpms = reverse sort {$a<=>$b} @tpms;

my $min_tpm_thresh = int($tpms[0]);
my $num_features = 1;

print "neg_min_tpm\tnum_features\n";

shift @tpms;
while (@tpms) {

    my $tpm = shift @tpms;

    if ($tpm < $min_tpm_thresh) {
        print "" . (-1*$min_tpm_thresh) . "\t$num_features\n";
        $min_tpm_thresh = int($tpm);

    }
    $num_features++;
}

print "$min_tpm_thresh\t$num_features\n";

exit(0);
