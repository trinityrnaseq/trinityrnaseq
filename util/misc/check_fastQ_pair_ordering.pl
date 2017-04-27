#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fastq_reader;

my $usage = "usage: $0 left.fq right.fq\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;

main: {

    my $left_fq_reader = new Fastq_reader($left_fq);
    my $right_fq_reader = new Fastq_reader($right_fq);

    my $ok_counter = 0;
    my $error_counter = 0;

    while (my $left_fq_record = $left_fq_reader->next()) {
        
        my $right_fq_record = $right_fq_reader->next() or die "Error, no next record in $right_fq";
        
        my $left_core_acc = $left_fq_record->get_core_read_name();
        my $left_full_acc = $left_fq_record->get_full_read_name();

        my $right_core_acc = $right_fq_record->get_core_read_name();
        my $right_full_acc = $right_fq_record->get_full_read_name();

        if ($left_core_acc eq $right_core_acc) {
            $ok_counter++;
        }
        else {
            $error_counter++;
            print join("\t", $left_full_acc, $right_full_acc, $left_core_acc, $right_core_acc, "ERROR") . "\n\n";
            print STDERR "\r[$ok_counter ok, $error_counter error]   ";
            

            die "Error, found out-of-order pairing. Stopping now.";
            
        }
        if ($ok_counter % 1000 == 0) {
            print STDERR "\r[$ok_counter ok, $error_counter error]   ";
        }
        
    }
    

    exit(0);
}
        
        
