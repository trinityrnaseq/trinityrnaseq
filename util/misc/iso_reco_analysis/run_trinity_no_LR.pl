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
    
    my $outdir = dirname($genome_file);
    
    my $cmd = "$ENV{TRINITY_HOME}/Trinity --seqType fa --single $outdir/simul.reads.fa --CPU 1 --max_memory 1G --output $outdir/trinity_no_LR_outdir --full_cleanup --trinity_complete";
    
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


