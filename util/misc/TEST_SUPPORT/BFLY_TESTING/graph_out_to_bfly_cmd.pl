#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 bfly_graph_out_list.file  bfly_jar\n\n";

my $bfly_graph_list_file = $ARGV[0] or die $usage;
my $bfly_jar = $ARGV[1] or die $usage;

main: {

    open (my $fh, $bfly_graph_list_file) or die $!;
    while (<$fh>) {
        chomp;
        my $graph_out_file = $_;

        $graph_out_file =~ s/\.out$//;

        unless (-s "$graph_out_file.reads") {
            print STDERR "WARNING - missing $graph_out_file.reads file.... skipping this one.\n";
            next;
        }
        
        my $cmd = "java -Xmx10G -Xms1G  -XX:ParallelGCThreads=2  -jar $bfly_jar -N 100000 -L 200 -F 500 -C $graph_out_file --path_reinforcement_distance=75 ";
                
        print "$cmd\n";
    }
    close $fh;


    exit(0);
}

