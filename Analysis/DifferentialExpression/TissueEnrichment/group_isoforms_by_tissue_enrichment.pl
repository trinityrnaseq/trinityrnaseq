#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 up_cats.counts\n\n";

my $counts_file = $ARGV[0] or die $usage;


main: {

    my %data;

    open (my $fh, "$counts_file") or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $trans = $x[0];
        my $category = $x[5];

        my $gene = $trans;
        $gene =~ s/\_i\d+$//;
        
        $data{$gene}->{$trans} = $category;
    }
    close $fh;


    foreach my $gene (keys %data) {

        my @trans = keys %{$data{$gene}};
        
        if (scalar @trans < 2) { next; }


        my %categories;
        
        foreach my $trans (@trans) {

            my $category = $data{$gene}->{$trans};

            $categories{$category}++;
        }

        my @cats = keys %categories;

        print join("\t", $gene, scalar(@trans), scalar @cats, join("__", sort @cats) ) . "\n";
    }


    exit(0);
}
        
