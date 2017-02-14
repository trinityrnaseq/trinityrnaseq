#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::Bin/../../PerlLib");
use Pipeliner;
use Cwd;


my $usage = "\n\n\tusage: $0 input.fasta reads.left.[fq|fa] reads.right.[fq|fa]\n\n";

my $target_fa = $ARGV[0] or die $usage;
my $left_reads = $ARGV[1] or die $usage;
my $right_reads = $ARGV[2] or die $usage;

unless ($target_fa =~ /^\//) {
    $target_fa = cwd() . "/$target_fa";
}

main: {

    my $pipeliner = new Pipeliner('-verbose' => 2);

    $pipeliner->add_commands(new Command("samtools faidx $target_fa", "$target_fa.fai.ok"));
    
    my $cmd = "bowtie2-build $target_fa $target_fa";
    $pipeliner->add_commands(new Command($cmd, "$target_fa.bowtie2-build.ok"));


    my $format = ($left_reads =~ /q$/i) ? '-q' : '-f';

    my $alignments_file = cwd() . "/alignments.$$.bam";
    
    $cmd = "set -eof pipefail; bowtie2 --no-unal -X 1000 -x $target_fa $format -1 $left_reads -2 $right_reads | samtools view -Sb - | samtools sort  > $alignments_file";
    $pipeliner->add_commands(new Command($cmd, "$alignments_file.ok"));

    $pipeliner->add_commands(new Command("samtools index $alignments_file", "$alignments_file.bai.ok"));

    $pipeliner->run();

    system("igv.sh -g $target_fa $alignments_file");


    exit(0);
}



