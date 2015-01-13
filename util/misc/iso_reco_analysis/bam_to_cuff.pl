#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "usage: $0 bams.list\n\n";

my $bam_list_file = $ARGV[0] or die $usage;

open (my $fh, $bam_list_file) or die $!;
while (<$fh>) {
    chomp;
    my $bam_file = $_;

    my $base_dir = dirname($bam_file);
    
    my $cmd = "cufflinks -o $base_dir/ $bam_file";

    print "$cmd\n";
}

exit(0);


