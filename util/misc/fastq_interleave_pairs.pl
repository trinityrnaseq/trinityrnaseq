#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use Fastq_reader;

my $usage = "\n\nusage: $0 left.fq right.fq\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;

main: {

    my $left_fastq_reader = new Fastq_reader($left_fq);
    my $right_fastq_reader = new Fastq_reader($right_fq);

    while (1) {
        my $left_fq = $left_fastq_reader->next();
        if ($left_fq) {
            print $left_fq->get_fastq_record();
        }


        my $right_fq = $right_fastq_reader->next();
        if ($right_fq) {
            print $right_fq->get_fastq_record();
        }

        if ($left_fq && $right_fq) {
            my $left_core_read_name = $left_fq->get_core_read_name();
            my $right_core_read_name = $right_fq->get_core_read_name();
            
            if ($left_core_read_name ne $right_core_read_name) {
                die "Error, read names are out of synch: $left_core_read_name vs. $right_core_read_name ";
            }
        }
        
        else {
            last;
        }
        
    }
    
    if ($left_fastq_reader->next() || $right_fastq_reader->next()) {
        die "Error, unequal number of fastq records in left and right fq files.";
    }


    exit(0);
}


        
