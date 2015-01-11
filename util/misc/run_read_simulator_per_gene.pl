#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use FindBin;


my $usage = "usage: $0 file.fasta [max_genes]\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $max_genes = $ARGV[1];


main: {

    my $sim_out_dir = "sim_AS_data";
    unless (-d $sim_out_dir) {
        mkdir $sim_out_dir or die $!;
    }

    my $fasta_reader = new Fasta_reader($fasta_file);

    my %gene_to_seqs;

    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        
        my ($trans, $gene) = split(/;/, $acc);
        
        unless ($gene) {
            die "Error, need trans;gene  format for accession: $acc";
        }

        my $sequence = $seq_obj->get_sequence();

        push (@{$gene_to_seqs{$gene}}, { acc => $acc,
                                         seq => $sequence, });
        
    }

   
    my $gene_counter = 0;
    ## only including those entries that are alt-spliced
    foreach my $gene (keys %gene_to_seqs) {

        my @trans = @{$gene_to_seqs{$gene}};

        if (scalar @trans == 1) {
            next;
        }
        

        my $outdir = $gene;
        $outdir =~ s/\W/_/g;
        
        mkdir ("$sim_out_dir/$outdir") or die $!;
        
        my $template_file = "$sim_out_dir/$outdir/$outdir.template.fa";
        open (my $ofh, ">$template_file") or die "Error, cannot write to $template_file";
        foreach my $entry (@trans) {
            my ($acc, $sequence) = ($entry->{acc}, $entry->{seq});
            print $ofh ">$acc\n$sequence\n";
        }
        close $ofh;
        
        my $outfile = "$sim_out_dir/$outdir/$outdir.reads.fa";
        
        my $cmd = "$FindBin::Bin/../simulate_illuminaPE_from_transcripts.pl --transcripts $template_file --SS > $outfile";
        &process_cmd($cmd);
        
        $gene_counter++;

        if ($max_genes && $gene_counter >= $max_genes) {
            last;
        }

        
    }

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


