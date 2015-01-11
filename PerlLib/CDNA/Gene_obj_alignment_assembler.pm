package main;
our $SEE;

=head1 NAME

CDNA::Gene_obj_alignment_assembler;

=cut

=head1 DESCRIPTION

This package assembles multiple gene objs by assembling compatible sets of coordinates.  This module inherits from the CDNA::Genome_based_cDNA_assembler module.  Please see this inherited module for additional methods available.

=cut

package CDNA::Gene_obj_alignment_assembler;
use strict;
use base qw (CDNA::PASA_alignment_assembler);
use CDNA::CDNA_alignment;


=item new()

=over 4

B<Description:> Instantiate a new CDNA::Gene_obj_alignment_assembler object.

B<Parameters:> $seq_ref

$seq_ref is a scalar reference to the genomic sequence string.

B<Returns:> obj_ref 

=back

=cut


sub new {
    my $packagename = shift;
    my $seqref = shift;
    my $self = $packagename->SUPER::new();
    $self->{sequence_ref} = $seqref;
    bless ($self, $packagename);
   
    return ($self);
}
    


=item assemble_genes()

=over 4

B<Description:> assembles a list of Gene_obj(s).  Each of the Gene_objs

B<Parameters:> @Gene_obj

@Gene_obj is an array of gene objects created via the Gene_obj.pm module.

B<Returns:> @CDNA_alignments

The @CDNA_alignments array contains the list of CDNA::CDNA_alignment objects created based on the gene objects.

The CDNA_alignment accession is set to the Model_feat_name of the gene_obj.  Retrieving the cDNA acc using the get_acc() method will allow a mapping back to the gene object.  Also, the newly created CDNA_alignment objects that are returned should be in the same order as the inputed gene objects, so indexing should provide mapping info as well.

use the get_assemblies() function to fetch the merged genes.

=back

=cut





sub assemble_genes {
    my $self = shift;
    my @gene_objs = @_;
    my @cDNA_alignments;
    foreach my $gene_obj (@gene_objs) {
	my $alignment_obj = $self->gene_obj_to_cdna_alignment($gene_obj);
	push (@cDNA_alignments, $alignment_obj);
    }
    $self->assemble_alignments(@cDNA_alignments);
    return (@cDNA_alignments);
}



=item gene_obj_to_cdna_alignment()

=over 4

B<Description:> method converts a Gene_obj to a CDNA::CDNA_alignment obj.

B<Parameters:> Gene_obj

B<Returns:> CDNA::CDNA_alignment

=back

=cut

sub gene_obj_to_cdna_alignment {
    my $self = shift;
    my $gene_obj = shift;
    my @exons = $gene_obj->get_exons();
    my %coords;
    my @alignment_segments;
    my $cdna_length = 0;
    foreach my $exon (@exons) {
	my ($end5, $end3) = $exon->get_coords();
	my $alignment_seg = new CDNA::Alignment_segment($end5, $end3);
	push (@alignment_segments, $alignment_seg);
	$cdna_length += abs ($end3 - $end5) + 1;
    }
    my $alignment_obj = new CDNA::CDNA_alignment($cdna_length, \@alignment_segments, $self->{sequence_ref});
    my $feat_name = $gene_obj->{Model_feat_name};
    print "Gene has the following model feat_name: $feat_name\n" if $SEE;
    $alignment_obj->set_acc($feat_name);
    $alignment_obj->set_title($gene_obj->{com_name});
    $alignment_obj->set_fli_status(1); #treat like a fli cdna
    my $align_orient = $alignment_obj->get_spliced_orientation();
    ## Make sure strand is set correctly (should only matter for single exon genes).
    if ($align_orient ne $gene_obj->{strand}) {
	$alignment_obj->set_orientation($gene_obj->{strand});
	$alignment_obj->set_spliced_orientation($gene_obj->{strand});
    }
    $alignment_obj->remap_cdna_segment_coords();
    print $alignment_obj->toToken() ." " . $alignment_obj->get_acc() .  "\n" if $SEE;
    return ($alignment_obj);
}


1; #EOM




