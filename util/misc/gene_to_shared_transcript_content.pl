#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use Nuc_translator;

my $usage = "usage:\t $0 refTrans.fasta [strandSpecific=0]\n\n";

my $trans_fa = $ARGV[0] or die $usage;

my $strand_specific_flag = $ARGV[1] or 0;

my $kmer_size = 24;

main: {


    my %kmer_to_gene;

    my $counter = 0;
    my $fasta_reader = new Fasta_reader($trans_fa);
    while (my $seq_obj = $fasta_reader->next()) {

        my $accession = $seq_obj->get_accession();
        my $sequence = uc $seq_obj->get_sequence();

        $counter++;
        print STDERR "\r[$counter] -tracking $accession";
        

        my ($trans, $gene) = split(/;/, $accession);
        unless ($gene) {
            die "Error, need format: >trans;gene  , instead found: $accession ";
        }
        
        for (my $i = 0; $i <= length($sequence) - $kmer_size; $i++) {
            
            my $kmer = substr($sequence, $i, $kmer_size);
            
            unless (length($kmer) == $kmer_size) {
                die "Error, didn't extract proper kmer length $kmer_size : [$kmer] ";
            }

            unless ($strand_specific_flag) {
                my @kmers = ($kmer, &reverse_complement($kmer));
                @kmers = sort @kmers;
                $kmer = $kmers[0];
            }

            $kmer_to_gene{$kmer}->{$gene} = 1;
        }
    }

    print STDERR "\n\n-reorganizing as gene sets to kmer lists.\n";
    
    my %genes_to_kmers;
    foreach my $kmer (keys %kmer_to_gene) {

        my @genes = keys %{$kmer_to_gene{$kmer}};

        if (scalar @genes > 1) {
            my $gene_token = join(",", sort @genes);
            push (@{$genes_to_kmers{$gene_token}}, $kmer);
        }
    }

    print STDERR "\n\n-outputting gene links and kmers.\n";
    
    foreach my $gene_set (sort keys %genes_to_kmers) {

        my @kmers = sort @{$genes_to_kmers{$gene_set}};

        print "$gene_set\t" . join(",", sort @kmers) . "\n";
    }


    exit(0);
}
