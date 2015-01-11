#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 (FPKM|RAW) fileA fileB ...\n\n";

unless (@ARGV && scalar(@ARGV) > 2 && $ARGV[0] =~ /^(FPKM|RAW)$/) {
    die $usage;
}

my $dat_type = shift @ARGV;
my @files = @ARGV;

my %data;

foreach my $file (@files) {

    open (my $fh, $file) or die "Error, cannot open file $file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { 
            next;
        }
        my @x = split(/\t/);
        my $trans = $x[0];
        
        my $val = ($dat_type eq "FPKM") ? $x[5] : $x[4];
        $data{$trans}->{$file} = $val;
    }
    close $fh;
}


print "#transcript\t" . join("\t", @files) . "\n";

foreach my $transcript (sort keys %data) {
    
    print $transcript;
    
    foreach my $file (@files) {
        
        my $val = $data{$transcript}->{$file} || 0;
        print "\t$val";
    }
    print "\n";
}

exit(0);

