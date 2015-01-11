#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 targets.fasta\n\n";

my $target_fasta_file = $ARGV[0] or die $usage;


main: {
    
    my $fasta_reader = new Fasta_reader($target_fasta_file);

    my %seqs = $fasta_reader->retrieve_all_seqs_hash();
    
    foreach my $acc (keys %seqs) {
	my $sequence = $seqs{$acc};

	my $seq_len = length($sequence);
	if ($seq_len < 500) { next; }

	my $bubble_missing_seq = $sequence;
	$bubble_missing_seq = substr($bubble_missing_seq, 0, 200) . substr($bubble_missing_seq, 350);
	
	my $new_gene_acc = $acc;
	$new_gene_acc =~ s/\W/_/g;

	print ">isoA-$new_gene_acc;$new_gene_acc\n$sequence\n"
	    . ">isoB-$new_gene_acc;$new_gene_acc\n$bubble_missing_seq\n";

    }
    

    exit(0);
}

