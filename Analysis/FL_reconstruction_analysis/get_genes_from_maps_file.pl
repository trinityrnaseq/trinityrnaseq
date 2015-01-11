#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.maps (G|T)\n\n";

my $maps_file = $ARGV[0] or die $usage;
my $GorT = $ARGV[1] or die $usage;


main: {

    open (my $fh, $maps_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($genes_info, $asms) = split(/\t/);

        my @genes = split(/,/, $genes_info);
        foreach my $gene_text (@genes) {
            my ($trans,$gene) = split(/;/, $gene_text);
            
            if ($GorT eq "G") {
                print "$gene\n";
            }
            else {
                print "$trans\n";
            }
        }
    }
    close $fh;

    exit(0);
}

