#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd;

my $usage = "usage: $0 genome_fa_files.list\n\n";

my $genome_fa_files = $ARGV[0] or die $usage;

my $workdir = cwd();


open (my $fh, $genome_fa_files) or die "Error, cannot open file $genome_fa_files";
while (<$fh>) {
    chomp;
    my $genome_file = $_;
    
    my $outdir = dirname($genome_file);
    
    if ($outdir !~ /^\./) {
        $outdir = "$workdir/$outdir";
    }
    

    my $cmd = "/home/unix/bhaas/GITHUB/trinityrnaseq/Trinity --seqType fa --single $outdir/simul.reads.fa --CPU 1 --max_memory 5G --output $outdir/trinity_WITH_LR_outdir --full_cleanup --long_reads $outdir/simul.transcriptome.cdnas";
    
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


