#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "usage: $0 genome_fa_files.list\n\n";

my $genome_fa_files = $ARGV[0] or die $usage;

open (my $fh, $genome_fa_files) or die "Error, cannot open file $genome_fa_files";
while (<$fh>) {
    chomp;
    my $genome_file = $_;
    
    my $gff3_file = $genome_file;
    $gff3_file =~ s/\.fa$/\.gff3/;

    my $outdir = dirname($gff3_file);

    my $cmd = "/home/unix/bhaas/GITHUB/trinityrnaseq/util/misc/simulate_reads_sam_and_fa.pl --gff3 $gff3_file --genome $genome_file --frag_length 300 --read_length 76 --SS_lib_type F --out_prefix $outdir/simul";
    
    #&process_cmd($cmd);
    print "$cmd\n";
    
}
exit(0);


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


