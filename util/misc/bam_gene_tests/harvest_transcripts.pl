#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "usage: $0 bfly_fasta_files.list.file\n\n";

my $bfly_list_file = $ARGV[0] or die $usage;


main: {

    open (my $fh, $bfly_list_file) or die $!;
    while (<$fh>) {
        chomp;
        
        my $file = $_;

        my $base = basename($file);
        my @pts = split(/\./, $base);
        my $core = shift @pts;

        open (my $fh2, $file) or die "Error, cannot open file $file";
        while (<$fh2>) {
            
            if (/>/) {
                s/>/>$core-/;
            }
            print;
        }
        close $fh2;
    }
    close $fh;

    exit(0);
}


            
