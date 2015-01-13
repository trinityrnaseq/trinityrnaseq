#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd;

my $usage = "usage: $0 trin_fa.list\n\n";

my $trin_fa_files = $ARGV[0] or die $usage;

my $workdir = cwd();


open (my $fh, $trin_fa_files) or die "Error, cannot open file $trin_fa_files";
while (<$fh>) {
    chomp;
    my $trin_fa_file = $_;
    
    my $outdir = dirname($trin_fa_file);

    my $cmd = "gmap -g $outdir/gene.fa $trin_fa_file -f 3 > $trin_fa_file.gff3";
    
    
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


