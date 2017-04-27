#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use List::Util qw(shuffle);


my $usage = <<__EOUSAGE__;

######################################################################################
#
# Required:
#
#  --cmds <string>                cmds file
#
#  --max_batch_size <int>         maximum batch size
#
# Optional:
#
#  --shuffle                      shuffle the commands in random order before batching
#
#######################################################################################


__EOUSAGE__

    ;


my $help_flag;

my $cmds_file;
my $max_batch_size;
my $shuffle_flag = 0;

&GetOptions ( 'h' => \$help_flag,
              'cmds=s' => \$cmds_file,
              'max_batch_size=i' => \$max_batch_size,
              'shuffle' => \$shuffle_flag);


if ($help_flag) {
    die $usage;
}

unless ($cmds_file && $max_batch_size) {
    die $usage;
}


main: {
    
    my @cmds = `cat $cmds_file`;
    chomp @cmds;

    if ($shuffle_flag) {
        @cmds = shuffle(@cmds);
    }

    my $num_cmds = scalar(@cmds);

    my $cmds_per_batch = int($num_cmds / $max_batch_size);
    if ($cmds_per_batch < 1) {
        $cmds_per_batch = 1;
    }
    
    while (@cmds) {
        
        my @batch;
        for (1..$cmds_per_batch) {
            my $cmd = shift @cmds;
            if ($cmd) { 
                push (@batch, $cmd);
            }
        }
        my $batched_cmds = join(" && ", @batch);
        
        print "$batched_cmds\n";
    }

    exit(0);
}



    
    


