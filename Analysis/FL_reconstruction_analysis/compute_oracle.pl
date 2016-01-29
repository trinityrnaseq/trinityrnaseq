#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $usage = "usage: $0 reads.fasta ref_transcripts.fasta [SS]\n\n";

my $reads_file = $ARGV[0] or die $usage;
my $ref_transcripts_fasta = $ARGV[1] or die $usage;
my $SS_flag = $ARGV[2] || 0;


my $cmd = "$FindBin::RealBin/../../Inchworm/bin/inchworm "
    . " --reads $reads_file "
    . " --checkFastaPath $ref_transcripts_fasta ";


unless ($SS_flag) {
    $cmd .= " --DS ";
}

my $ret = system($cmd);
if ($ret) {
    die "Error, CMD: $cmd died with ret $ret";
}

exit(0);

