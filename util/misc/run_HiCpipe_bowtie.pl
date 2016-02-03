#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


my $usage = "usage: $0 genome.fasta  left.fq  right.fq [output_dir]\n\n";

my $genome_file = $ARGV[0] or die $usage;
my $left_fq_file = $ARGV[1] or die $usage;
my $right_fq_file = $ARGV[2] or die $usage;
my $output_dir = $ARGV[3] || "bowtie.$$.dir";

while ($output_dir =~ m|/$|) {
    chop $output_dir;
}


main: {

    ## run bowtie
    my $cmd = "$FindBin::RealBin/../alignReads.pl --target $genome_file --left $left_fq_file --right $right_fq_file "
        . " --seqType fq --aligner bowtie -o $output_dir  --max_dist_between_pairs 900000000 --no_rsem --retain_intermediate_files "
        . " -- -a -m 1 --best --strata -p 4 --chunkmbs 512 ";
    &process_cmd($cmd) unless (-s "$output_dir/$output_dir.nameSorted.sam");

    $cmd = "$FindBin::RealBin/HiCpipe_nameSortedSam_to_raw.pl $output_dir/$output_dir.nameSorted.sam > $output_dir/$output_dir.raw";
    &process_cmd($cmd) unless (-s "$output_dir/$output_dir.raw");
    

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
