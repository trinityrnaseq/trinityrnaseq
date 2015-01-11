#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 fileA fileB ...\n\n";

my @files = @ARGV;
unless (@files) {
    die $usage;
}

main: {

    my %data;
    foreach my $file (@files) {
        open (my $fh, $file) or die $!;
        while (<$fh>) {
            chomp;
            unless (/\w/) { next; }
            my ($key, $rest) = split(/\t/, $_, 2);
            $rest =~ s/\t/ /g;
            $data{$key}->{$file} = $rest;
        }
        close $fh;
    }

    print join("\t", @files) . "\n";
    foreach my $acc (keys %data) {
        print "$acc";
        foreach my $file (@files) {
            my $val = $data{$acc}->{$file};
            unless (defined $val) {
                $val = "NA";
            }
            print "\t$val";
        }
        print "\n";
    }

    exit(0);

}
        
