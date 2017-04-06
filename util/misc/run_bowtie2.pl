#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Carp;


my $CPU = 2;

my $usage = <<__EOUSAGE__;

############################################################################
#
#  --target <string>        target for alignment
#
#  --left <string>          read_1.fq
#
#  optional:
#
#  --right <string>         read_2.fq
#
#  --CPU <int>              number of threads (default: $CPU)
#
#

 usage: $0  --target target.seq --left reads_1.fq [--right reads_2.fq --CPU 8]
    
      and you can pipe it into samtools to make a bam file:
   
         | samtools view -Sb - | samtools sort - -o bowtie2.coordSorted.bam

#############################################################################



__EOUSAGE__

    ;

    
my $help_flag;
my $target_seq;
my $reads_1_fq;
my $reads_2_fq;


&GetOptions ( 'h' => \$help_flag,
              'target=s' => \$target_seq,
              'left=s' => \$reads_1_fq,
              'right=s' => \$reads_2_fq,
              'CPU=i' => \$CPU,
    );


if ($help_flag) { die $usage; }

unless ($target_seq && $reads_1_fq) { die $usage; }



main: {

    unless (-s "$target_seq.1.bt2") {
        my $cmd = "bowtie2-build --threads $CPU $target_seq $target_seq 1>&2 ";
        &process_cmd($cmd);
    }
    
    my $format = ($reads_1_fq =~ /\.fq|\.fastq/) ? "-q" : "-f";
    
    my $bowtie2_cmd = "bowtie2 --local --no-unal -k 50 -x $target_seq $format ";
    if ($reads_2_fq) {
        $bowtie2_cmd .= " -1 $reads_1_fq -2 $reads_2_fq ";
    }
    else {
        $bowtie2_cmd .= " -U $reads_1_fq ";
    }

    
    &process_cmd($bowtie2_cmd);
    
    exit(0);
}
    
