#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 cmds.list  cmds.list.cache_success\n\n";


my $cmds_list_file = $ARGV[0] or die $usage;
my $cached_successes_file = $ARGV[1] or die $usage;



main: {


    my %OK;
    
    {
        # get cached successes
        
        open (my $fh, $cached_successes_file) or die $!;
        while (<$fh>) {
            chomp;
            my $cmd = $_;

            $OK{$cmd} = 1;
        }
        close $fh;

    }

    open (my $fh, $cmds_list_file) or die $!;
    while (<$fh>) {
        chomp;
        my $cmd = $_;
        unless ($OK{$cmd}) {
            print "$cmd\n";
        }
    }

    close $fh;

    exit(0);
}



