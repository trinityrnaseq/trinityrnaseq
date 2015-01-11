#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.maps\n";

my $file = $ARGV[0] or die $usage;

my %genes;
my %isoforms;

my %gene_merges;
my %isoform_merges;


open (my $fh, $file) or die $!;


while (<$fh>) {
	chomp;
	my @x = split(/\t/);
	my $genes = shift @x;
	
	
	my @entries = split(/,/, $genes);
	
	
	my %iso;
	my %g;
	
	foreach my $entry (@entries) {
		
		$iso{$entry}++;
		
		$entry =~ s/^.*;//;
	
		$g{$entry}++;
	}

	my $num_genes = scalar (keys %g);

	
	my ($gene_counter_href, $isoform_counter_href) = ($num_genes > 1) 
		? (\%gene_merges, \%isoform_merges) 
		: (\%genes, \%isoforms);
	
	foreach my $gene (keys %g) {
		$gene_counter_href->{$gene}++;
	}
	foreach my $isoform (keys %iso) {
		$isoform_counter_href->{$isoform}++;
	}
	
	
}


my $num_single_genes = scalar(keys %genes);
my $num_single_isoforms = scalar(keys %isoforms);

my $num_merged_genes = scalar(keys %gene_merges);
my $num_isoform_merges = scalar(keys %isoform_merges);


#print "#file\tsingle_gene\tsingle_isoforms\tmerged_genes\tmerged_isoforms\n";
print "$file\t$num_single_genes\t$num_single_isoforms\t$num_merged_genes\t$num_isoform_merges\n";


exit(0);

