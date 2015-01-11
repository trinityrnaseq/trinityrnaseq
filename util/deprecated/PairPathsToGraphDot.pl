#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib", "$FindBin::Bin/../PerlLib/KmerGraphLib");

use StringGraph;
use StringNode;
use ColorGradient;
use Data::Dumper;
use List::Util qw(shuffle);

my $usage = "usage: $0 Butterfly.verbose.log\n\n";

my $butterfly_log = $ARGV[0] or die $usage;


main: {

    my @pair_paths;
    
    {
        open (my $fh, $butterfly_log) or die $!;
        while (<$fh>) {
            
            if (/PAIRPATH: PairPath \[_paths=\[\[([^\]]+)\], \[([^\]]*)\]\]\]=(\d+)/) {
                my $path_pt1 = $1;
                my $path_pt2 = $2;
                my $count = $3;
                push (@pair_paths, [$path_pt1, $path_pt2, $count]);
               
            }
            else {
                print "IGNORING: $_";
            }
        }
        close $fh;

    }

    unless (@pair_paths) {
        die "Error, no pair paths found.";
    }

    my $graph = new StringGraph();
    
    my @colors = &ColorGradient::get_RGB_gradient(scalar(@pair_paths));
    @colors = &ColorGradient::convert_RGB_hex(@colors);
    
    @colors = shuffle(@colors);

    my $counter = 0;

    foreach my $pair_path (@pair_paths) {
        
        print STDERR Dumper($pair_path);
        
        $counter++;

        my $color = shift @colors;
        
        my ($path_pt1, $path_pt2, $count) = @$pair_path;
        
        my @path_node_names = &get_path_node_names($path_pt1);

        print STDERR "Path 1: @path_node_names\n";
        
        if (scalar(@path_node_names) == 1) {
            @path_node_names = (@path_node_names, @path_node_names); # include a self edge so evidence will show.
        }
        

        $graph->add_sequence_to_graph($counter, \@path_node_names, $counter, $color);
        
        if ($path_pt2) {
            my @path_node_names_2 = &get_path_node_names($path_pt2);
    
            print STDERR "Path 2: @path_node_names_2\n";
        
            if (scalar(@path_node_names_2) == 0) {
                @path_node_names_2 = (@path_node_names_2, @path_node_names_2);
            }
            
            $graph->add_sequence_to_graph($counter, \@path_node_names_2, $counter, $color);

        }
    }
    
    print $graph->toGraphViz();
    
}

####
sub get_path_node_names {
    my ($path_text) = @_;

    
    my @node_names = split(/[\s,]+/, $path_text);
    
    return(@node_names);
}

