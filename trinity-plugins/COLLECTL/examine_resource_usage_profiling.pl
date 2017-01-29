#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;

my $usage = "usage: $0 collectl_dir/ \n\n";

my $collectl_dir = $ARGV[0] or die $usage;


main: {

    my $cmd = "$FindBin::Bin/util/collectl_dat_to_time_matrix.py --dat $collectl_dir/collectl.dat";
    &process_cmd($cmd);

    $cmd = "$FindBin::Bin/util/plot_time_vs_resource.Rscript $collectl_dir/collectl";
    &process_cmd($cmd);

    exit(0);
}


    
