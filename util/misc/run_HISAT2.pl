#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = "\n\n\tusage: $0 genome.fa genome.gtf samples_file.txt\n\n";

my $genome_fa = $ARGV[0] or die $usage;
my $annotation_gtf = $ARGV[1] or die $usage;
my $samples_file = $ARGV[2] or die $usage;

main: {

    my @read_sets = &parse_samples_file($samples_file);
    
    ## align reads to the mini-genome using hisat2

    unless (-e "$genome_fa.hisat2.build.ok") {
    
        &process_cmd("hisat2_extract_splice_sites.py $annotation_gtf  > $annotation_gtf.ss");
        
        &process_cmd("hisat2_extract_exons.py $annotation_gtf > $annotation_gtf.exons");
        
        &process_cmd("hisat2-build --exon $annotation_gtf.exons --ss $annotation_gtf.ss  $genome_fa $genome_fa");

        &process_cmd("touch $genome_fa.hisat2.build.ok");
    }
    
    foreach my $read_set_aref (@read_sets) {
        my ($sample_id, $left_fq, $right_fq) = @$read_set_aref;

        
        &process_cmd("set -eof pipefail; hisat2 --dta -x $genome_fa -1 $left_fq -2 $right_fq | samtools view -Sb -F 4 | samtools sort -o $sample_id.cSorted.bam");

    }
    
    exit(0);
    
}



####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples;
    
    open(my $fh, $samples_file) or die "Error, cannot open file $samples_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($cond, $rep, $fq_a, $fq_b) = @x;

        push (@samples, [$rep, $fq_a, $fq_b]);
    }
    close $fh;

    return (@samples);
}

