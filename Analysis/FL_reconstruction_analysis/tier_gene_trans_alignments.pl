#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

my $usage = "usage: $0 genes.fasta transcripts.fasta max_intron BLAT_CPU [STRAND-SPECIFIC_FLAG=0]\n\n";

my $genes_fasta = $ARGV[0] or die $usage;
my $trans_fasta = $ARGV[1] or die $usage;
my $max_intron = $ARGV[2] or die $usage;
my $blat_cpu = $ARGV[3] or die $usage;
my $SS = $ARGV[4] || 0;


main: {

    ## run blat:
    my $cmd = "$FindBin::RealBin/../../util/process_BLAT_alignments.pl -g $genes_fasta -t $trans_fasta -I $max_intron --CPU $blat_cpu --KEEP_PSLX";

    &process_cmd($cmd);

    $cmd = "cat blat_out_dir/*top_1 > blat.top_1.pslx";
    &process_cmd($cmd);

    $cmd = "$FindBin::RealBin/util/blat_top_tier_genes.pl blat.top_1.pslx $SS";
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
