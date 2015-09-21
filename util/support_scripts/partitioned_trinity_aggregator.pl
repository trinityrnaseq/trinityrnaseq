#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 token < trinity_fasta_files_listing\n\n";

my $token = $ARGV[0] or die $usage;

while (<STDIN>) {
    my $filename = $_;
    chomp $filename;
    unless (-e $filename) {
        print STDERR "ERROR, filename: $filename is indicated to not exist.\n";
        next;
    }
    if (-s $filename) {
        
        # ie. read_partitions/Fb_0/CBin_161/c16167.trinity.reads.fa
        
        $filename =~ m|c(\d+)\.trinity\.reads| or die "Error, cannot parse Trinity component value from filename: $filename";
        my $component = $1;
        
        open (my $fh, $filename) or die "Error, cannot open file $filename";
        while (<$fh>) {
            if (/>/) {
                s/>/>${token}${component}\_/;
            }
            print;
        }
    }
}

exit(0);

