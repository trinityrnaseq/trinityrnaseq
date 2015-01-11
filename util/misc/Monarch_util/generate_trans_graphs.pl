#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 transcripts.cdna.fasta\n\n";

my $transcripts_fasta = $ARGV[0] or die $usage;


my $graphs_per_dir = 100;


main: {
    
    my $fasta_reader = new Fasta_reader($transcripts_fasta) or die $!;
    
    
    my $count = 0;
   
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $accession = $seq_obj->get_accession();
        $accession =~ s/\W/_/g;

        my $sequence = $seq_obj->get_sequence();

        my $dir_no = int($count/$graphs_per_dir);
        my $outdir = "trans_graphs/g_$dir_no";
        if (! -d $outdir) {
            &process_cmd("mkdir -p $outdir");
        }
        
        my $fa_file = "$outdir/$accession.fa";
        open (my $ofh, ">$fa_file") or die $!;
        print $ofh ">$accession\n$sequence\n";
        close $ofh;

        my $cmd = "~/SVN/trinityrnaseq/trunk/util/misc/Monarch --misc_seqs $fa_file --graph $fa_file.dot";
        &process_cmd($cmd);

        
        $count++;
    }


    print STDERR "\n\nDone.\n\n";
    
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
