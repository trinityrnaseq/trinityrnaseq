#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");

use Fastq_reader;

my $usage = "\n\nusage: $0 interleaved.fq  [left_output_filename.fq right_output_filename.fq]\n\n";

my $interleaved_fq = $ARGV[0] or die $usage;
my $left_out_filename = $ARGV[1] || "unweaved.left.$$.fq";
my $right_out_filename = $ARGV[2] || "unweaved.right.$$.fq";

main: {

    my $fastq_reader = new Fastq_reader($interleaved_fq);

    open (my $left_ofh, ">$left_out_filename") or die "Error, cannot write to $left_out_filename";
    open (my $right_ofh, ">$right_out_filename") or die "Error, cannot write to $right_out_filename";
    
    while (my $fq = $fastq_reader->next()) {

        my $read_name = $fq->get_full_read_name();
        
        my $record = $fq->get_fastq_record();

        if ($read_name =~ m|/1$|) {
            print $left_ofh $record;
        }
        elsif ($read_name =~ m|/2$|) {
            print $right_ofh $record;
        }
        else {
            die "Error, cannot decipher left or right read from name: $read_name";
        }
    }

    close $left_ofh;
    close $right_ofh;


    print "Done.\n";
    

    exit(0);
}


        
