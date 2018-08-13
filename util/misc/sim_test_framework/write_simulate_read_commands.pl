#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use File::Basename;

my $usage = "\n\n\tusage: $0 target_trans_files.list [opts ex. --wgsim ...]\n\n";

my $target_trans_files_file = $ARGV[0] or die $usage;
shift @ARGV;


main: {

    my @ref_files = `cat $target_trans_files_file`;
    chomp @ref_files;

    foreach my $file (@ref_files) {
        my $outdir = dirname($file);
        my $cmd = "$FindBin::Bin/run_simulate_reads.wgsim.pl -R $file -O $outdir @ARGV";
                
        print "$cmd\n";
    }

    exit(0);
}
    
