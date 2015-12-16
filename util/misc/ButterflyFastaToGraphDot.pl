#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib", "$FindBin::RealBin/../../PerlLib/KmerGraphLib");

use Fasta_reader;
use StringGraph;
use StringNode;
use ColorGradient;

my $usage = "usage: $0 Butterfly.fasta\n\n";

my $butterfly_fasta_file = $ARGV[0] or die $usage;


main: {

    my @seqs_n_paths;
    
    {
        open (my $fh, $butterfly_fasta_file) or die $!;
        while (<$fh>) {
            if (/^>/) {
                /^>(\S+) .* path=\[(.*)\]/ or die "Error, cannot parse header: $_";
                
                my $acc = $1;
                my $path = $2;
                push (@seqs_n_paths, [$acc, $path]);
            }
        }
        close $fh;

    }

    my $graph = new StringGraph();
    
    my @colors = &ColorGradient::get_RGB_gradient(scalar(@seqs_n_paths));
    @colors = &ColorGradient::convert_RGB_hex(@colors);
    
    foreach my $seq_n_path (@seqs_n_paths) {
        my ($acc, $path_text) = @$seq_n_path;
        
        my @path_node_names = &get_path_node_names($path_text);
        
        my $color = shift @colors;

        $graph->add_sequence_to_graph($acc, \@path_node_names, 1, $color);
        
    }
    
    print $graph->toGraphViz();
    
}

####
sub get_path_node_names {
    my ($path_text) = @_;

    my @node_names;
    
    my @parts = split(/\s+/, $path_text);
    foreach my $part (@parts) {
        my ($node_name, $coords) = split(/:/, $part);
        
        my ($lend, $rend) = split(/-/, $coords);
        my $length = $rend-$lend + 1;
        
        $node_name = "${node_name}_len$length";
        
        push (@node_names, $node_name);
    }
    
    return(@node_names);
}

