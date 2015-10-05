#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 token < trinity_fasta_files_listing\n\n";

my $token = $ARGV[0] or die $usage;


my $counter = 0;

while (<STDIN>) {
    my $filename = $_;
    chomp $filename;
    unless (-e $filename) {
        print STDERR "ERROR, filename: $filename is indicated to not exist.\n";
        next;
    }
    if (-s $filename) {
        $counter++;
        open (my $fh, $filename) or die "Error, cannot open file $filename";
        while (<$fh>) {
            if (/>/) {
                s/>/>${token}_${counter}\_/;
            }
            print;
        }
    }
}

exit(0);

