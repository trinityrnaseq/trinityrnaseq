#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 sampleA.EXPRESS.isoform.results sampleB.EXPRESS.isoform.results ...\n\n";

unless (@ARGV) {
    die $usage;
}

my @eXpress_files = @ARGV;

if (scalar @eXpress_files == 1) {

    if (-s $eXpress_files[0]) {
        # allow for a file listing the various files.
        @eXpress_files = `cat $eXpress_files[0]`;
        chomp @eXpress_files;
    }
    else {
        die $usage;
    }
}


=header_format

0       bundle_id
1       target_id
2       length
3       eff_length
4       tot_counts
5       uniq_counts
6       est_counts
7       eff_counts
8       ambig_distr_alpha
9       ambig_distr_beta
10      fpkm
11      fpkm_conf_low
12      fpkm_conf_high
13      solvable

=cut


main: {
    
    my %data;
    
    foreach my $file (@eXpress_files) {
        print STDERR "-reading file: $file\n";
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>; # ignore it
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[1];
            my $count = $x[7];
            $data{$acc}->{$file} = $count;
        }
        close $fh;
    }
    
    my @filenames = @eXpress_files;
    foreach my $file (@filenames) {
        $file = basename($file);
    }

    print STDERR "\n\n* Outputting combined matrix.\n\n";
    
    
    print join("\t", "", @filenames) . "\n";
    foreach my $acc (keys %data) {
        
        print "$acc";

        foreach my $file (@eXpress_files) {

            my $count = $data{$acc}->{$file};
            unless (defined $count) {
                $count = "NA";
            }

            print "\t$count";
            
        }
        
        print "\n";
        
    }
    

    print STDERR "Done.\n\n";
    
    exit(0);
}
    
        
