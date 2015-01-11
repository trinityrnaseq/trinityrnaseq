package main;
our $SEE;

package AlignGraph;

use strict;
use warnings;
use Carp;
use AlignNode;

use base qw (ReadCoverageGraph);

no warnings qw (recursion);


sub new {
	my $packagename = shift;
	
	my $self = $packagename->SUPER::new();

	
	## Nodes exist in order of end5 -> end3 at all times.  Strand is included in the node name just for safety reasons.

	bless ($self, $packagename);

	return($self);
}


sub add_alignment {
	my $self = shift;
	
	my ($read_acc, $scaffold, $strand, $genome_coords_aref) = @_;

	my @align_positions;
	
	foreach my $coordset (sort {$a->[0]<=>$b->[0]} @$genome_coords_aref) {

		my ($lend, $rend) = sort {$a<=>$b} @$coordset;
		
		for (my $i = $lend; $i <= $rend; $i++) {
			push (@align_positions, $i);
		}
	}
	
	$self->_add_ordered_positions_to_graph(\@align_positions, $strand, $read_acc);
	
	return;

}



sub _add_ordered_positions_to_graph {
	my $self = shift;
	my ($ordered_positions_aref, $strand, $read_acc) = @_;

	unless (ref $ordered_positions_aref eq 'ARRAY') {
		confess "Error, require ordered position list";
	}
	unless ($strand =~ /^[\+\-]$/) {
		confess "strand must be: + or - ";
	}
	
	my @align_positions = @$ordered_positions_aref;

	if ($strand eq '-') {
		@align_positions = reverse @align_positions;
	}

	my $prev_align_pos = shift @align_positions;
	while (@align_positions) {
		my $next_align_pos = shift @align_positions;
		
		# print "\tadding $next_align_pos\n";
		
		my $prev_align_node = $self->get_or_create_node("$prev_align_pos,$strand", $read_acc);
		my $next_align_node = $self->get_or_create_node("$next_align_pos,$strand", $read_acc);
		
		$self->link_adjacent_nodes($prev_align_node, $next_align_node);
		
		$prev_align_pos = $next_align_pos;

	}
	
	return;
}

sub get_all_nodes {
	my $self = shift;

	return( sort {$a->{_value} cmp $b->{_value}} $self->SUPER::get_all_nodes());

}



1; #EOM


=CIGAR_format

from: http://bioperl.org/pipermail/bioperl-l/2003-March/011591.html

cigar line format (where CIGAR stands for Concise
Idiosyncratic Gapped Alignment Report).

In the cigar line format alignments are sotred as follows:

M: Match
D: Deletino
I: Insertion

An example of an alignment for a hypthetical protein match is shown 
below:


Query:   42 PGPAGLP----GSVGLQGPRGLRGPLP-GPLGPPL...
             PG    P    G     GP   R      PLGP
Sbjct: 1672 PGTP*TPLVPLGPWVPLGPSSPR--LPSGPLGPTD...


protein_align_feature table as the following cigar line:

7M4D12M2I2MD7M



 From SAM documentation:

Clipped alignment. In Smith-Waterman alignment, a sequence may not be aligned from the first residue to the last one.
Subsequences at the ends may be clipped off. We introduce operation ʻSʼ to describe (softly) clipped alignment. Here is
an example. Suppose the clipped alignment is:

REF: AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATCGATGTG
READ: gggGTGTAACC-GACTAGgggg

where on the read sequence, bases in uppercase are matches and bases in lowercase are clipped off. The CIGAR for
this alignment is: 3S8M1D6M4S.
Spliced alignment. In cDNA-to-genome alignment, we may want to distinguish introns from deletions in exons. We
introduce operation ʻNʼ to represent long skip on the reference sequence. Suppose the spliced alignment is:

REF: AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATCGATGTG
READ: GTGTAACCC................................TCAGAATA

where ʻ...ʼ on the read sequence indicates the intron. The CIGAR for this alignment is: 9M32N8M.

=cut

