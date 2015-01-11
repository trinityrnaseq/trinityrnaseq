#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use SAM_reader;
use SAM_entry;

use Fasta_reader;
use CIGAR;

my $usage = "usage: $0 file.sam genome.fa min_len=10\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $genome_file = $ARGV[1] or die $usage;
my $min_len = $ARGV[2] || 10;

main: {

	my $genome_sref;
	my $prev_scaff = "";

	my $fasta_reader = new Fasta_reader($genome_file);
	my %genome = $fasta_reader->retrieve_all_seqs_hash();
	
	my $sam_reader = new SAM_reader($sam_file);
	
	while ($sam_reader->has_next()) {
		
		my $sam_entry = $sam_reader->get_next();
		
		my $cigar = $sam_entry->get_cigar_alignment();
		
		my $sequence = $sam_entry->get_sequence();
		
		my $scaff = $sam_entry->get_scaffold_name();
		
		my $strand = $sam_entry->get_query_transcribed_strand();
		
		if ($scaff ne $prev_scaff) {
		    my $contig = $genome{$scaff};
			$genome_sref = \$contig;
		}
		
		my ($genome_align_aref, $query_align_aref) = $sam_entry->get_alignment_coords();

				
		my @genome_coords = @$genome_align_aref;
		my @query_coords = @$query_align_aref;

		if (&get_len($genome_coords[0]) < $min_len) {
			shift @genome_coords;
			shift @query_coords;
		}
		if (&get_len($genome_coords[$#genome_coords]) < $min_len) {
		    pop @genome_coords;
			pop @query_coords;
		}
		
		my $new_cigar = &CIGAR::construct_cigar(\@genome_coords, \@query_coords, length($sequence), 
												$genome_sref, $strand);
		
		#print "$cigar\t$new_cigar\n" if ($cigar ne $new_cigar);
	
		my $new_genome_start = $genome_coords[0]->[0];
		
		my @fields = $sam_entry->get_fields();
		$fields[3] = $new_genome_start;
		$fields[5] = $new_cigar;
		
		print join("\t", @fields) . "\n";
		
	}
	

	exit(0);
}


####
sub get_len {
	my ($seg) = @_;

	my ($lend, $rend) = @$seg;
	
	my $len = $rend - $lend + 1;

	return($len);
}
