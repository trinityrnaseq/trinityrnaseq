#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;


my $usage = "usage: $0 transcripts.cdna.fasta\n\n";

my $transcripts_fasta_file = $ARGV[0] or die $usage;

my $graphs_per_dir = 100;


main: {
    
    my $fasta_reader = new Fasta_reader($transcripts_fasta_file) or die $!;
    
    
    my $count = 0;
   
    my %gene_to_seq;

    while (my $seq_obj = $fasta_reader->next()) {
        
        my $accession = $seq_obj->get_accession();

        my $sequence = $seq_obj->get_sequence();
        
        my ($trans, $gene) = split(/;/, $accession);


        $gene_to_seq{$gene}->{$trans} = $sequence;

    }

    foreach my $gene (keys %gene_to_seq) {
        
        my $trans_href = $gene_to_seq{$gene};
        
        my @trans = keys %$trans_href;
        unless (scalar @trans > 1) {
            # want just alt-splice ones.
            next; 
        }

        $gene =~ s/\W/_/g;
        
        
        my $dir_no = int($count/$graphs_per_dir);
        my $outdir = "gene_altSplice_graphs/g_$dir_no";
        if (! -d $outdir) {
            &process_cmd("mkdir -p $outdir");
        }

        my $fa_file = "$outdir/$gene.fa";
        open (my $ofh, ">$fa_file") or die $!;
        foreach my $trans_acc (keys %$trans_href) {

            my $seq = $trans_href->{$trans_acc};
            print $ofh ">$trans_acc\n$seq\n";
            
        }
        
        
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
