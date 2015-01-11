#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $usage = "usage: $0 bfly.A.fasta [bfly.B.fasta ...]\n\n";

my @files = @ARGV or die $usage;

main: {

    my %data;
    
    foreach my $file (@files) {
        open (my $fh, $file) or die "Error, cannot open file $file";
        while (<$fh>) {
            
            unless (/^>/) { next; }
            
            my $comp_id;
            if (/>\S*(c\d+\.graph_c\d+)_seq\d+/) {
                $comp_id = $1;
            }
            elsif (/^>\S*(comp\d+_c\d+)_seq\d+/) {
                $comp_id = $1;
            }
            else {
                die "Error, couldn't extract component identifier from $_";
            }
        
            $data{$comp_id}->{$file}++;
        }
    }

    ## output data
    print join("\t", "#component", @files) . "\n";
    foreach my $component (keys %data) {
        
        print "$component";
        foreach my $file (@files) {
            my $count = $data{$component}->{$file} || 0;
            print "\t$count";
        }
        print "\n";
    }

    exit(0);
}


            

