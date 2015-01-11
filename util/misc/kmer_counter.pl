#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Ktree;
use Nuc_translator;

my $usage = "usage: $0 fasta_file KmerSize [DSmode]\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $kmer_length = $ARGV[1] or die $usage;
my $DS_mode = $ARGV[2] || 0;


main: {
    
    my $ktree = new Ktree();

    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $seq = $seq_obj->get_sequence();
        my @chars = split(//, $seq);

        for (my $i = 0; $i <= length($seq) - $kmer_length; $i++) {
            
            my $kmer = join("", @chars[$i..($i+$kmer_length-1)]);
            
            $ktree->add_kmer($kmer);
            
            if ($DS_mode) {
                $ktree->add_kmer( &reverse_complement($kmer) );
            }
        }
    }
    
    $ktree->report_kmer_counts();
    
    exit(0);
}


            
