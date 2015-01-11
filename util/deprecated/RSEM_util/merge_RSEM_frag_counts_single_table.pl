#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 sampleA.RSEM.isoform.results sampleB.RSEM.isoform.results ...\n\n";

unless (@ARGV) {
    die $usage;
}

my @rsem_files = @ARGV;

if (scalar @rsem_files == 1) {

    if (-s $rsem_files[0]) {
        # allow for a file listing the various files.
        @rsem_files = `cat $rsem_files[0]`;
        chomp @rsem_files;
    }
    else {
        die $usage;
    }
}


=header_format

0       transcript_id
1       gene_id
2       length
3       effective_length
4       expected_count
5       TPM
6       FPKM
7       IsoPct

=cut


main: {


    my %data;
    
    foreach my $file (@rsem_files) {
        print STDERR "-reading file: $file\n";
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>; # ignore it
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[0];
            my $count = $x[4];
            $data{$acc}->{$file} = $count;
        }
        close $fh;
    }
    
    my @filenames = @rsem_files;
    foreach my $file (@filenames) {
        $file = basename($file);
    }

    print STDERR "\n\n* Outputting combined matrix.\n\n";
    
    
    print join("\t", "", @filenames) . "\n";
    foreach my $acc (keys %data) {
        
        print "$acc";

        foreach my $file (@rsem_files) {

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
    
        
