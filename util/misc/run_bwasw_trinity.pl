#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "usage: $0 Trinity.fasta reads.{fa,fq}\n\n";

my $trinity_fasta = $ARGV[0] or die $usage;
my $reads_fasta = $ARGV[1] or die $usage;


main: {

    my $cmd = "fasta_file_header_stripper.pl < $trinity_fasta > $trinity_fasta.noheader";
    &process_cmd($cmd);

    $cmd = "bwa index -a is $trinity_fasta.noheader";
    &process_cmd($cmd);

    $cmd = "samtools faidx $trinity_fasta.noheader";
    &process_cmd($cmd);

    my $outfile_prefix = basename($reads_fasta);
    $cmd = "bwa bwasw $trinity_fasta.noheader $reads_fasta > $outfile_prefix.sam";
    &process_cmd($cmd);

    $cmd = "samtools view -bt $trinity_fasta.noheader.fai $outfile_prefix.sam > $outfile_prefix.bam";
    &process_cmd($cmd);

    $cmd = "samtools sort $outfile_prefix.bam -o $outfile_prefix.bam ";
    &process_cmd($cmd);

    $cmd = "samtools index $outfile_prefix.bam";
    &process_cmd($cmd);

    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}
