#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 file.kmers bfly.log\n\n";

my $file_kmers = $ARGV[0] or die $usage;
my $bfly_log = $ARGV[1] or die $usage;

my %kmers;
{
    open (my $fh, $file_kmers) or die $!;
    while (<$fh>) {
        chomp;
        $kmers{$_} = 1;
    }
    close $fh;
}

open (my $fh, $bfly_log) or die "Error, cannot open file $bfly_log";
while (<$fh>) {
    my $line = $_;
    chomp;
    # EDGE_PRUNING::removeLightOutEdges() removing the edge: G:W-1(V6691_D-1) GAGGCTGTGAAGAGACTGGCAGAG -> G:W-1(V9713_D-1) GAGGCTGTGAAGAGACTGGCAGAG (weight: 1.0 <= e_edge_thr: 6.550000000000001, EDGE_THR=0.05

    if (/^EDGE_PRUNING/) {
        my @x = split(/\s+/);
        my $kmer_A = $x[5];
        my $kmer_B = $x[6];

        if ($kmers{$kmer_A} && $kmers{$kmer_B}) {
            print "!!\t$line";
        }
    }

}
exit(0);


            
    
