#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 DE.graph\n\n";

my $DE_graph = $ARGV[0] or die $usage;


=color_panel

> colorpanel(10, 'black', 'purple', 'red')
 [1] "#000000" "#28083C" "#501078" "#7818B4" "#A020F0" "#A020F0" "#B818B4"
 [8] "#D01078" "#E7083C" "#FF0000"

=cut


my @colors = (
    "#000000", "#28083C", "#501078", "#7818B4", "#A020F0", "#A020F0", "#B818B4",
    "#D01078", "#E7083C", "#FF0000"
    );


main: {

    print "digraph G {\n";
    
    open (my $fh, $DE_graph) or die "Error, cannot open file $DE_graph";
    my $line = <$fh>;
    chomp $line;
    close $fh;

    my @x = split(/\t/, $line);
    shift @x;
    foreach my $pair (@x) {
        my ($from, $to, $logFC) = split(/,/, $pair);
        
        my $color_index = int($logFC + 0.5);
        if ($color_index > $#colors) {
            $color_index = $#colors;
        }
        my $color = $colors[$color_index];
        
        print "    $from->$to\[color=\"$color\"]\n";
    }
    
    print "}\n";
    
    exit(0);
}


        
