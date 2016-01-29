#!/usr/bin/env perl

use strict;
use warnings;


use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use Fasta_reader;

my $usage = "usage: $0 fasta_file\n\n";

my $fasta_file = $ARGV[0] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);


while (my $seq_obj = $fasta_reader->next()) {

	my $acc = $seq_obj->get_accession();

	my $contig_acc = $acc;

	$acc =~ s/;/_/;

	my $seq = $seq_obj->get_sequence();

	my $seq_len = length($seq);


	print join("\t", $contig_acc, ".", "transcript", 1, $seq_len, ".", "+", ".", 
			   "gene_id \"g.$acc\"; transcript_id \"t.$acc\";") . "\n";

	print join("\t", $contig_acc, ".", "exon", 1, $seq_len, ".", "+", ".", 
			   "gene_id \"g.$acc\"; transcript_id \"t.$acc\";") . "\n";

	print "\n";
}


exit(0);


