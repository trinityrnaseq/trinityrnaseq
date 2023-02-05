#!/usr/bin/env perl

package main;
our $DEBUG;

package Gene_obj;
use strict;
use Nuc_translator;
#use Gene_ontology;
use Longest_orf;
use Storable qw (store retrieve freeze thaw dclone);
use warnings;
use Data::Dumper;
use Carp qw (croak cluck confess);
#use URI::Escape;

=head1 NAME

package Gene_obj

=cut



=head1 DESCRIPTION

    Gene_obj(s) encapsulate the elements of both gene structure and gene function. The gene structure is stored in a hierarchical fashion as follows:

    Gene  =========================================================

    Exon  =========     =========         =========        ========

    CDS      ======     =========         ======

   
    where a Gene is a container for Exon(s), and each Exon is a container for a CDS, and an Exon can contain a single CDS component.  An Exon lacking a CDS exon is an untranslated exon or UTR exon.  The region of an Exon which extends beyond the CDS is also considered a UTR.
  
    
    There are several ways to instantiate gene objects.  A simple example is described:

    Exon and CDS component coordinates can be assigned as hashes.

    ie. 
    
    my %mrna = ( 100 => 200,
	         300 => 500 );

    my %CDS = ( 150=>200,
		300=>450);

    my $sequence = "GACTACATTTAATAGGGCCC"; #string representing the genomic sequence
    my $gene = new Gene_obj();
    
    $gene->{com_name} = "hypothetical protein";

    $gene->populate_gene_obj(\%CDS, \%mRNA, \$sequence);
    print $gene->toString();

    
    
    Alternatively, the individual components of genes (Exons and CDSs) can be instantiated separately and used to build the Gene from the ground up (See packages mRNA_exon_obj and CDS_exon_obj following this Gene_obj documentation).
    
    my $cds_exon = new CDS_exon_obj (150, 200);
    
    my $mRNA_exon = new mRNA_exon_obj (100, 200);
   
    $mRNA_exon->set_CDS_exon_obj($cds_exon);

    my $gene_obj = new Gene_obj ();

    $gene_obj->{gene_name} = "hypothetical gene";
    $gene_obj->{com_name} = "hypothetical protein";
  
    $gene_obj->add_mRNA_exon_obj($mRNA_exon);

    $gene_obj->refine_gene_object();

    $gene_obj->create_all_sequence_types (\$sequence);  #ref to genomic sequence string.    

    print $gene_obj->toString();


    The API below describes useful functions for navigating and manipulating the Gene object along with all of its attributes.
    


=cut






=over 4

=item new()

B<Description:> Constructor for Gene_obj 

B<Parameters:> none

B<Returns:> $gene_obj


The Gene_obj contains several attributes which can be manipulated directly (or by get/set methods if they exist).  These attributes include:

    asmbl_id # identifier for the genomic contig for which this gene is anchored.
    TU_feat_name #feat_names are TIGR temporary identifiers.
    Model_feat_name # temp TIGR identifier for gene models
    locus  #identifier for a gene (TU) ie. T2P3.5
    pub_locus  #another identifier for a gene (TU)   ie. At2g00010
    model_pub_locus #identifier for a gene model (model)  ie. At2g00010.1
    model_locus #analagous to locus, but for model rather than gene (TU)
    alt_locus   #alternative locus 
    gene_name # name for gene
    com_name  # name for gene product 
    comment #internal comment
    pub_comment #comment related to gene
    ec_num   # enzyme commission number
    gene_sym  # gene symbol
    is_5prime_partial # 0|1  missing start codon.
    is_3prime_partial # 0|1  missing stop codon.
    is_pseudogene # 0|1
    curated_com_name # 0|1
    curated_gene_structure # 0|1
    
    ## Other attributes set internally  Access-only, do not set directly.
        
    gene_length  # length of gene span (int).
    mid_pt  # holds midpoint of gene-span
    strand  # [+-]
    protein_seq # holds protein sequence
    protein_seq_length
    CDS_sequence  #holds CDS sequence (translated to protein); based on CDS_exon coordinates
    CDS_seq_length 
    cDNA_sequence  #holds cDNA sequence; based on mRNA exon coordinates.
    cDNA_seq_length 
    gene_sequence #holds unspliced transcript
    gene_sequence_length #length of unspliced transcript
    gene_type # "protein-coding", #default type for gene object.  Could be changed to "rRNA|snoRNA|snRNA|tRNA" to accommodate other gene or feature types.
    num_additional_isoforms # int 
    
    
=back

=cut



sub new {
    shift;
    my $self = { asmbl_id => 0, #genomic contig ID
                 locus => undef,       #text
                 pub_locus => undef,   #text  ie. At2g00010
                 model_pub_locus =>undef, #text ie. At2g00010.1
                 model_locus => undef, #text ie. F12G15.1
                 alt_locus => undef,   #text
                 gene_name => undef, #text
                 com_name => undef,    #text
                 comment => undef,
                 curated_com_name => 0,
                 curated_gene_structure => 0,
                 pub_comment => undef, #text
                 ec_num => undef, #text (enzyme commission number)
                 gene_type => "protein-coding", #default type for gene object.  Could be changed to "rRNA|snoRNA|snRNA|tRNA" to accomodate other gene or feature types.
                 gene_sym => undef, #text (gene symbol)
                 mRNA_coords => 0, #assigned to anonymous hash of end5->end3 relative to the parent sequence
                 CDS_coords => 0,  #assigned to anonymous hash of end5->end3 relative to the parent sequence
                 mRNA_exon_objs => 0,  # holds arrayref to mRNA_obj, retrieve only thru method: get_exons()
                 num_exons => 0,      # number of exons in this gene_obj
                 model_span => [],     # holds array ref to (end5,end3) for CDS range of gene.
                 gene_span => [],      # holds array ref to (end5,end3) for mRNA range of gene.
                 gene_length => 0,     # length of gene span (int).
                 mid_pt => 0,         # holds midpoint of gene-span
                 strand => 0,      # [+-]
                 gi => undef,          #text
                 prot_acc => undef,     #text
                 is_pseudogene => 0, # toggle indicating pseudogene if 1.
                 is_5prime_partial => 0, #boolean indicating missing 5' part of gene.
                 is_3prime_partial => 0, #boolean
                 protein_seq => undef,    # holds protein sequence
                 protein_seq_length => 0,
                 CDS_sequence => undef,    #holds CDS sequence (translated to protein); based on CDS_exon coordinates
                 CDS_seq_length => 0,
                 cDNA_sequence => undef,   #holds cDNA sequence; based on mRNA exon coordinates.
                 cDNA_seq_length => 0,
                 gene_sequence => undef, #holds unspliced transcript
                 gene_sequence_length => 0, #length of unspliced transcript
                 TU_feat_name => undef,    #feat_names are TIGR temporary identifiers.
                 Model_feat_name =>undef,
                 classification => 'annotated_genes', #type of seq_element.
                 gene_synonyms => [],    #list of synonymous model feat_names
                 GeneOntology=>[], #list of Gene_ontology assignment objects.  ...see GeneOntology.pm
                 
                 ## Additional functional attributes:
                 secondary_gene_names => [],
                 secondary_product_names => [],
                 secondary_gene_symbols => [],
                 secondary_ec_numbers =>[],
                 
                 
                 ## Alternative splicing support.  
                 num_additional_isoforms => 0,  # number of additional isoforms stored in additonal_isoform list below
                 additional_isoforms => [] # stores list of Gene_objs corresponding to the additional isoforms.
		     
                 };
    bless($self);
    return ($self);
}




=over 4

=item erase_gene_structure()

B<Description:> Removes the structural components of a gene (ie. exons, CDSs, coordinate spans, any corresponding sequences)

B<Parameters:> none

B<Returns:> none 

=back

=cut


## erase gene structure
sub erase_gene_structure {
    my $self = shift;
    $self->{mRNA_exon_objs} = 0;
    $self->{num_exons} = 0;
    $self->{model_span} = [];
    $self->{gene_span} = [];
    $self->{gene_length} = 0;
    $self->{strand} = 0;
    $self->{protein_seq} = 0;
    $self->{CDS_sequence} = 0;
    $self->{CDS_seq_length} = 0;
    $self->{cDNA_sequence} = 0;
    $self->{cDNA_seq_length} = 0;
}


=over 4

=item clone_gene()

B<Description:> Clones this Gene_obj by copying attributes from this Gene to a new gene.  Does NOT do a deep clone for all attributes.  See dclone() for a more rigorous cloning method.  This method is safer because all references are not cloned, only the critical ones.

B<Parameters:> none

B<Returns:> new Gene_obj

=back

=cut



## all objects are cloned.  References to data only are not.
sub clone_gene {
    my $self = shift;
    my $clone = new Gene_obj();
    
    
    ## Copy over the non-ref attribute values.
    foreach my $key (keys %$self) {
        my $value = $self->{$key};
        if (defined $value) {
            ## Not copying over refs.
            if (ref $value) {
                next;
            }
            
            ## Not copying over attributes of length > 200, such as protein/nucleotide sequences
            my $length = length($value);
            if ($length > 200) { next;}
        }
        
        # passed tests above, copying attribute.
        $clone->{$key} = $value;
        
    }
    
    ## copy over the gene synonyms.
    my @gene_syns = @{$self->{gene_synonyms}};
    $clone->{gene_synonyms} = \@gene_syns;
    
    
    ## copy the GO assignments:
    my @GO_assignments = $self->get_gene_ontology_objs();
    if (@GO_assignments) {
        foreach my $go_assignment (@GO_assignments) {
            my $go_clone = dclone($go_assignment);
            $clone->add_gene_ontology_objs($go_clone);
        }
    }
    
    
    ## copy gene structure.
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        $clone->add_mRNA_exon_obj($exon->clone_exon());
    }
    
    foreach my $isoform ($self->get_additional_isoforms()) {
        my $isoform_clone = $isoform->clone_gene();
        $clone->add_isoform($isoform_clone);
    }
    
    $clone->refine_gene_object();
    
    return ($clone);
}




=over 4

=item deep_clone()

B<Description:> Provides a deep clone of a gene_obj.  Only references supported in Gene_obj documentation are supported.  Those added in a rogue way are undef()d

B<Parameters:> none

B<Returns:> $gene_obj

uses the Storable dclone() function to deep clone the Gene_obj

=back

=cut


    ;
## all objects are cloned.  References to data only are not.
sub deep_clone {
    my $self = shift;
    my $clone = dclone($self);
    
    my %supported_refs = (model_span => 1,
                          gene_span => 1,
                          gene_synonyms => 1,
                          Gene_ontology => 1,
                          additional_isoforms=>1,
                          mRNA_exon_objs => 1);
    
    foreach my $gene_obj ($clone, $clone->get_additional_isoforms()) {
        
        my @keys = keys %$gene_obj;
        foreach my $key (@keys) {
            my $value = $gene_obj->{$key};
            if (ref $value && !$supported_refs{$key}) {
                $gene_obj->{$key} = undef;
            }
        }
    }
    
    return ($clone);
}


=over 4

=item populate_gene_obj()

B<Description:> Given CDS and mRNA coordinates stored in hash form, a gene object is populated with mRNA and CDS exons.  This is one available way to populate a newly instantiated Gene_obj.

B<Parameters:> $cds_hash_ref, $mRNA_hash_ref, <$seq_ref>

$mRNA_hash_ref is a reference to a hash holding the end5 => end3 coordinates of the Exons

$cds_hash_ref same as mRNA_has_ref except holds the CDS end5 => end3 coordinates.

$seq_ref is a reference to a string containing the genomic sequence.  This is an optional parameter.


B<Returns:> none

=back

=cut

    ;

## Do several things at once: assign CDS and mRNA coordinates, and build gene sequences.
## The \$seq_ref is optional in case you want to create the sequence types.
sub populate_gene_obj {
    my ($self, $cds_ref, $mRNA_ref, $seq_ref) = @_;
    $self->set_CDS_coords ($cds_ref);
    $self->set_mRNA_coords ($mRNA_ref);
    $self->refine_gene_object();
    if (ref $seq_ref) {
        $self->create_all_sequence_types($seq_ref);
    }
    ## reinitialize the hashrefs:
    $self->{mRNA_coords} = 0;
    $self->{CDS_coords} = 0;
    
    
}


# alias above
sub populate_gene_object {
    my $self = shift;
    $self->populate_gene_obj(@_);
}


####
sub populate_gene_object_via_CDS_coords {
	my $self = shift;
	my @coordsets = @_;

	foreach my $coordset (@coordsets) {
		my ($end5, $end3) = @$coordset;
		my $mrna_exon_obj = mRNA_exon_obj->new($end5, $end3);
		my $cds_obj = CDS_exon_obj->new($end5, $end3);
		$mrna_exon_obj->{CDS_exon_obj} = $cds_obj;
		$self->add_mRNA_exon_obj($mrna_exon_obj);
	}
	
	$self->refine_gene_object();
	return;
}


sub build_gene_obj_exons_n_cds_range {
	my $self = shift;
	my ($exons_aref, $cds_lend, $cds_rend, $orient) = @_;

	my @exon_coords;
	foreach my $exon_aref (@$exons_aref) {
		my ($exon_lend, $exon_rend) = sort {$a<=>$b} @$exon_aref;
		push (@exon_coords, [$exon_lend, $exon_rend] );
	}
	@exon_coords = sort {$a->[0]<=>$b->[0]} @exon_coords;


	unless ($orient =~ /^[\+\-]$/) {
		confess "Error, orient not [+-] ";
	}

	## build the CDS coordinates.

	my @cds_range;

    if ($cds_lend > 0 && $cds_rend > 0) {
        
        ($cds_lend, $cds_rend) = sort {$a<=>$b} ($cds_lend, $cds_rend);
        
        foreach my $exon_coords_aref (@exon_coords) {
            my ($exon_lend, $exon_rend) = @$exon_coords_aref;
            
            if ($exon_rend >= $cds_lend && $exon_lend <= $cds_rend) {
                
                ## got overlap
                my $cds_exon_lend = ($cds_lend < $exon_lend) ? $exon_lend : $cds_lend;
                
                my $cds_exon_rend = ($cds_rend > $exon_rend) ? $exon_rend : $cds_rend;
                
                push (@cds_range, [$cds_exon_lend, $cds_exon_rend]);
            }
        }
        
        unless (@cds_range) {
            confess "Error, no CDS exon coords built based on exon overlap";
        }
    }
	## all coordinate sets are ordered left to right.
	# build the coordinates href
	
	my %exon_coords;
	my %cds_coords;
	foreach my $exon_coords_aref (@exon_coords) {
		my ($exon_lend, $exon_rend) = @$exon_coords_aref;
		my ($exon_end5, $exon_end3) = ($orient eq '+') ? ($exon_lend, $exon_rend) : ($exon_rend, $exon_lend);
		$exon_coords{$exon_end5} = $exon_end3;
	}
	foreach my $cds_coords_aref (@cds_range) {
		my ($cds_lend, $cds_rend) = @$cds_coords_aref;
		my ($cds_end5, $cds_end3) = ($orient eq '+') ? ($cds_lend, $cds_rend) : ($cds_rend, $cds_lend);
		$cds_coords{$cds_end5} = $cds_end3;
	}

	# print Dumper (\%cds_coords) . Dumper (\%exon_coords);
	
	$self->populate_gene_obj(\%cds_coords, \%exon_coords);

	return ($self);
}


####
sub join_adjacent_exons {
    my $self = shift;

    my @exons = $self->get_exons();
    
    my $strand = $self->get_orientation();
    
    my $first_exon = shift @exons;
    my @new_exons = ($first_exon);

    while (@exons) {
        my $prev_exon = $new_exons[$#new_exons];
        my ($prev_end5, $prev_end3) = $prev_exon->get_coords();

        my $next_exon = shift @exons;
        my ($next_end5, $next_end3) = $next_exon->get_coords();

        if ( ($strand eq '+' && $prev_end3 == $next_end5 - 1)  # adjacent
             ||
             ($strand eq '-' && $prev_end3 == $next_end5 + 1) ) {
            
            $prev_exon->merge_exon($next_exon);
        }
        else {
            push (@new_exons, $next_exon);
        }
    }
    
    $self->{mRNA_exon_objs} = [@new_exons];
    
    $self->refine_gene_object();

    return;
}




=over 4

=item AAToNucleotideCoords()

B<Description:> Converts an amino acid -based coordinate to a genomic sequence -based coordinate.

B<Parameters:> $aa_coord

B<Returns:> $genomic_coord

undef is returned if the aa_coord could not be converted.


=back

=cut

    ;

sub AAToNucleotideCoords{
    my($self) = shift;
    my($aacoord) = shift;
    my($debug) = shift;
    my($PCDS_coords) = {};
    my($A2NMapping) = {};
    my($currAA) = 1; 
    my $strand = $self->{strand};
    my @exons = $self->get_exons();
    my($cds_count)=0;
    my($translated_bp)=-1;
    my($lastcarryover)=0; 
    my($end_bp);
    foreach my $exon (sort {
        if($strand eq "+"){
            $a->{end5}<=>$b->{end5};
        }
        else{
            $b->{end5}<=>$a->{end5};
        }
    } @exons) {
        my $cds = $exon->get_CDS_obj();
        if ($cds) {
            my @cds_coords = $cds->get_CDS_end5_end3();
            my($bpspread) = abs($cds_coords[0]-$cds_coords[1]);
            $bpspread+=$lastcarryover;
            my($nextAA) = int($bpspread/3); # last complete AA in CDS
            $lastcarryover = $bpspread%3;
            $PCDS_coords->{$currAA} = $currAA+$nextAA-1;
            if($strand eq "+"){
                $A2NMapping->{$currAA} = $cds_coords[0]<$cds_coords[1]?$cds_coords[0]:$cds_coords[1];
            }
            else{
                $A2NMapping->{$currAA} = $cds_coords[0]<$cds_coords[1]?$cds_coords[1]:$cds_coords[0];
            }
            print "DEBUG: $strand $cds_count AA range ($currAA - $PCDS_coords->{$currAA}) nucleotide start($A2NMapping->{$currAA})\n" if($debug);
            $currAA = $currAA+$nextAA;
            $cds_count++;
            if($strand eq "+"){
                $end_bp = $cds_coords[0]<$cds_coords[1]?$cds_coords[1]:$cds_coords[0];
            }
            else{
                $end_bp = $cds_coords[0]<$cds_coords[1]?$cds_coords[0]:$cds_coords[1];
            }
        }
    }
    # PCDS_coords key/value are start/stop aa counts for each cds;
    # A2NMapping stores cds AA start key to cds nucleotide start
    $cds_count=0;
    foreach my $PCDS_end5 (sort {
        $a<=>$b;
	}(keys %$PCDS_coords)) {
        my($PCDS_end3) = $PCDS_coords->{$PCDS_end5};
	    if($aacoord>=$PCDS_end5 && $aacoord<=$PCDS_end3){
            my($nucleotide_start) = $A2NMapping->{$PCDS_end5}; 
            my($aa_offset) = $aacoord - $PCDS_end5;
            my($nucleotide_offset) = $aa_offset*3;
            print "DEBUG: CDS offset $aa_offset AA $nucleotide_offset bp\n" if($debug);
            if($strand eq "+"){
                $translated_bp = $nucleotide_start+$nucleotide_offset;
            }
            else{
                $translated_bp = $nucleotide_start-$nucleotide_offset;
            }
            print "DEBUG: Mapping $aacoord to $translated_bp in cds $cds_count\n" if($debug);
            print "DEBUG: CDS $PCDS_end5 - $PCDS_end3 nucleotide start $A2NMapping->{$PCDS_end5}, nuc offset $nucleotide_offset\n" if($debug); 
	    }
        
        $cds_count++;
	}
    #}
    if($translated_bp == -1){
        $translated_bp = undef;
        print STDERR "Unable to translate AA coordinate: $aacoord. Off end. Using undef\n" if($debug);
    }
    return $translated_bp;
}



## private method, used by populate_gene_obj()
# sets CDS_coords instance member to a hash reference of CDS coordinates.   $hash{end5} = end3
sub set_CDS_coords {
    my $self = shift;
    my $hash_ref = shift;
    if (ref ($hash_ref) eq 'HASH') {
        $self->{CDS_coords} = $hash_ref;
    } else {
        print STDERR "Cannot set CDS_coords, must have hash reference\n";
    }
}




=over 4

=item get_gene_span()

B<Description:> Retrieves the coordinates which span the length of the gene along the genomic sequence.

B<Parameters:> none

B<Returns:> (end5, end3)

These coordinates represent the minimal and maximal exonic coordinates of the gene.  Orientation can be inferred by the relative values of end5 and end3.


=back

=cut

    ;

## All return gene end5, end3 ###
sub get_gene_span {
    my $self = shift;
    return (@{$self->{gene_span}});
}




## private
sub get_seq_span {
    my $self = shift;
    return ($self->get_gene_span());
}



=over 4

=item get_coords()

B<Description:> See get_gene_span()

B<Parameters:> none

B<Returns:> (end5, end3)

=back

=cut


sub get_coords {
    my $self = shift;
    return ($self->get_gene_span());
}



=over 4

=item get_model_span()

B<Description:> Retrieves the coordinates spanned by the protein-coding region of the gene along the genomic sequence.

B<Parameters:> none

B<Returns:> (end5, end3)

These coordinates are determined by the min and max of the CDS components of the gene.

=back

=cut




sub get_model_span {
    my $self = shift;
    return (@{$self->{model_span}});
}  


sub get_CDS_span { # preferred
	my $self = shift;
	return($self->get_model_span());
}


=over 4

=item get_transcript_span()

B<Description:> Retrieves the coordinates spanned by the exonic regions of the gene along the genomic sequence.

B<Parameters:> none

B<Returns:> (lend, rend)

These coordinates are determined by the min and max of the CDS components of the gene.

=back

=cut


sub get_transcript_span {
	my $self = shift;
	
	my @coords;
	my @exons = $self->get_exons();
	foreach my $exon (@exons) {
		push (@coords, $exon->get_coords());
	}
	@coords = sort {$a<=>$b} @coords;

	my $lend = shift @coords;
	my $rend = pop @coords;

	return($lend, $rend);
}


sub is_pseudogene {
    my $self = shift;
    return ($self->{is_pseudogene});
}

sub set_pseudogene {
    my $self = shift;
    my $pseudogene_val = shift;
    unless ($pseudogene_val =~ /[01]/) {
        confess "Error, can set pseudogene to zero or one only.\n";
    }

    foreach my $gene ($self, $self->get_additional_isoforms()) {
        $gene->{is_pseudogene} = $pseudogene_val;
    }

    return;
}



#private
# sets mRNA_coords instance member to a hash reference of CDS coordinates.   $hash{end5} = end3
sub set_mRNA_coords {
    my $self = shift;
    my $hash_ref = shift;
    if (ref ($hash_ref) eq 'HASH') {
        $self->{mRNA_coords} = $hash_ref;
    } else {
        print STDERR "Cannot set CDS_coords, must have hash reference\n";
    }
}


=over 4

=item refine_gene_object()

B<Description:> This method performs some data management operations and should be called at any time modifications have been made to the gene structure (ie. exons added or modified, model isoforms added, etc).  It performs the following orientations:

    -Sets (or resets)  gene span and model span coordinates, strand orientation, gene length, mid-point.

B<Parameters:> none

B<Returns:> none

=back

=cut

## Once mRNA_coords and CDS_coords have been assigned, this will populate the remaining elements in the gene object.

sub refine_gene_object {
    my ($self) = shift;
    #check to see if mRNA_coords field is populated.  If not, initialize.
    if ($self->{mRNA_coords} == 0) {
        $self->{mRNA_coords} = {};
    }
    my ($CDS_coords, $mRNA_coords) = ($self->{CDS_coords},  $self->{mRNA_coords});
    
    unless ($CDS_coords && $mRNA_coords) {
        #maybe created exon objects already
        if ($self->{mRNA_exon_objs}) {
            $self->trivial_refinement();
        }
        return;
    }
    # intialize mRNA_exon_objs to array ref.
    $self->{mRNA_exon_objs} = [];
    #retrieve coordinate data.
    my %mRNA = %$mRNA_coords;
    my %CDS = %$CDS_coords;
    my @mRNAcoords = keys %mRNA;
    my @CDScoords = keys %CDS;
    my (%new_mRNA, %new_CDS);
    ## if correlation between mRNA exons and CDS exons, then map CDS's to mRNA's, otherwise, replicate CDSs as mRNAs
    if ($#mRNAcoords >= $#CDScoords) {
        
        foreach my $mRNA_end5 (keys %mRNA) {
            my $mRNA_end3 = $mRNA{$mRNA_end5};
            #find overlapping cds exon to mRNA exon
            #easy to compare if in same orientation for all comparisons
            my ($m1, $m2) = ($mRNA_end5 < $mRNA_end3) ? ($mRNA_end5, $mRNA_end3) : ($mRNA_end3, $mRNA_end5);
            #create mRNA_exon_obj
            my $mRNA_exon_obj = mRNA_exon_obj->new ($mRNA_end5, $mRNA_end3);
            $new_mRNA{$mRNA_end5} = $mRNA_end3;
            foreach my $CDS_end5 (keys %CDS) {
                my $CDS_end3 = $CDS{$CDS_end5};
                my ($c1, $c2) = ($CDS_end5 < $CDS_end3) ? ($CDS_end5, $CDS_end3) : ($CDS_end3, $CDS_end5);
                ## do overlap comparison; CDS must be contained within mRNA exon
                if ( ($c1 >= $m1) && ($c2 <= $m2)) {
                    # found the contained CDS
                    $mRNA_exon_obj->{CDS_exon_obj} = CDS_exon_obj->new ($CDS_end5, $CDS_end3); 
                    $new_CDS{$CDS_end5} = $CDS_end3;
                    last;
                }
            }
            $self->add_mRNA_exon_obj($mRNA_exon_obj);
        }
    } else { # remap CDSs to mRNAS
        print STDERR "ERROR: mRNA exons < CDS exons.  Copying all CDS exons into mRNA exons. \n\n";
        foreach my $CDS_end5 (keys %CDS) {
            my $CDS_end3 = $CDS{$CDS_end5};
            my $mRNA_exon_obj = mRNA_exon_obj->new ($CDS_end5, $CDS_end3);
            $mRNA_exon_obj->{CDS_exon_obj} = CDS_exon_obj->new ($CDS_end5, $CDS_end3); 
            $self->add_mRNA_exon_obj($mRNA_exon_obj);
            $new_mRNA{$CDS_end5} = $CDS_end3;
            $new_CDS{$CDS_end5} = $CDS_end3;
        }
    } 
    
    $self->trivial_refinement();
    
    ## assign orientation to all children exon and CDS components.
    my $strand = $self->get_orientation();
    foreach my $exon ($self->get_exons()) {
        $exon->{strand} = $strand;
        if (my $cds = $exon->get_CDS_exon_obj()) {
            $cds->{strand} = $strand;
        }
    }
    return;
    
}


## alias
sub refine_gene_obj {
    my $self = shift;
    $self->refine_gene_object();
}

    
=over 4

=item get_exons()

B<Description:>Retrieves a list of exons belonging to this Gene_obj 

B<Parameters:> none

B<Returns:> @exons

@exons is an ordered list of mRNA_exon_obj; the first exon of the list corresponds to the first exon of the spliced gene.

=back

=cut

    ;

sub get_exons {
    my ($self) = shift;
    if ($self->{mRNA_exon_objs} != 0) {
        my @exons = (@{$self->{mRNA_exon_objs}});
        @exons = sort {$a->{end5}<=>$b->{end5}} @exons;
        if ($self->{strand} eq '-') {
            @exons = reverse (@exons);
        }
        return (@exons);
    } else {
        my @x = ();
        return (@x); #empty array 
    }
}


## private
sub get_segments {
    my $self = shift;
    return ($self->get_exons());
}



=over 4

=item number_of_exons()

B<Description:> Provides the number of exons contained by the Gene

B<Parameters:> none

B<Returns:> int

=back

=cut



sub number_of_exons {
    my $self = shift;
    my $exon_number = $#{$self->{mRNA_exon_objs}} + 1;
    return ($exon_number);
}







=over 4

=item get_intron_coordinates()

B<Description:> Provides an ordered list of intron coordinates

B<Parameters:> none

B<Returns:> ( [end5,end3], ....) 

A list of arrayRefs are returned providing the coordinates of introns, ordered from first intron to last intron within the gene.

=back

=cut

    ;

sub get_intron_coordinates {
    my $gene_obj = shift;
    my $strand = $gene_obj->get_orientation();
    my @exons = $gene_obj->get_exons();
    ## exon list should already be sorted.
    my @introns = ();
    
    my $num_exons = $#exons + 1;
    if ($num_exons > 1) { #only genes with multiple exons will have introns.
        if ($strand eq '+') {
            my $first_exon = shift @exons;
            while (@exons) {
                my $next_exon = shift @exons;
                my ($first_end5, $first_end3) = $first_exon->get_coords();
                my ($next_end5, $next_end3) = $next_exon->get_coords();
                my $intron_end5 = $first_end3 + 1;
                my $intron_end3 = $next_end5 -1;
                if ($intron_end5 < $intron_end3) {
                    push (@introns, [$intron_end5, $intron_end3]);
                }
                $first_exon = $next_exon;
            }
        } elsif ($strand eq '-') {
            my $first_exon = shift @exons;
            while (@exons) {
                my $next_exon = shift @exons;
                my ($first_end5, $first_end3) = $first_exon->get_coords();
                my ($next_end5, $next_end3) = $next_exon->get_coords();
                my $intron_end5 = $first_end3 - 1;
                my $intron_end3 = $next_end5 +1;
                if ($intron_end5 > $intron_end3) {
                    push (@introns, [$intron_end5, $intron_end3]);
                }
                $first_exon = $next_exon;
            }
            
        } else {
            die "Strand for gene_obj is not specified." . $gene_obj->toString();
        }
    }
    return (@introns);
}





#private
sub trivial_refinement {
    my $self = shift;
    my @exons = $self->get_exons();
    $self->{num_exons} = scalar(@exons);
    my (%mRNAexons, %CDSexons);
    foreach my $exon (@exons) {
        my ($exon_end5, $exon_end3) = $exon->get_mRNA_exon_end5_end3();
        $mRNAexons{$exon_end5} = $exon_end3;
        my $cds;
        if ($cds = $exon->get_CDS_obj()) {
            my ($cds_end5, $cds_end3) = $cds->get_CDS_end5_end3();
            $CDSexons{$cds_end5} = $cds_end3;
        }
    }
    my @mRNAexonsEnd5s = sort {$a<=>$b} keys %mRNAexons;
    my @CDSexonsEnd5s = sort {$a<=>$b} keys %CDSexons;
    my $strand = 0; #initialize.
    foreach my $mRNAend5 (@mRNAexonsEnd5s) {
        my $mRNAend3 = $mRNAexons{$mRNAend5};
        if ($mRNAend5 == $mRNAend3) {next;}
        $strand = ($mRNAend5 < $mRNAend3) ? '+':'-';
        last;
    }
    $self->{strand} = $strand;
    
    ## determine gene and model boundaries:
    my ($gene_end5, $gene_end3, $model_end5, $model_end3);
    my @gene_coords = sort {$a<=>$b} %mRNAexons;
    my @model_coords = sort {$a<=>$b} %CDSexons;
    my $gene_lend = shift @gene_coords;
    my $gene_rend = pop @gene_coords;
    ## bound gene by transcript span
    ($gene_end5, $gene_end3) = ($strand eq "+") ? ($gene_lend, $gene_rend) : ($gene_rend, $gene_lend);
    if (@model_coords) {
        ## bound model by protein coding span
        my $model_lend = shift @model_coords;
        my $model_rend = pop @model_coords;
        ($model_end5, $model_end3) = ($strand eq "+") ? ($model_lend, $model_rend) : ($model_rend, $model_lend);
    } else {
        ## give it gene boundaries instead:
        ($model_end5, $model_end3) = ($gene_end5, $gene_end3);
    }
    
    $self->{gene_span} = [$gene_end5, $gene_end3];
    $self->{gene_length} = abs ($gene_end3 - $gene_end5) + 1;
    $self->{mid_pt} = int (($gene_end5 + $gene_end3)/2);
    $self->{model_span} = [$model_end5, $model_end3]; 
    
    ## Refine isoforms if they exist.
    if (my @isoforms = $self->get_additional_isoforms()) {
        my @gene_span_coords = $self->get_gene_span();
        foreach my $isoform (@isoforms) {
            $isoform->refine_gene_object();
            push (@gene_span_coords, $isoform->get_gene_span());
        }
        @gene_span_coords = sort {$a<=>$b} @gene_span_coords;
        my $lend = shift @gene_span_coords;
        my $rend = pop @gene_span_coords;
        my $strand = $self->{strand};
        if ($strand eq '-') {
            ($lend, $rend) = ($rend, $lend);
        }
        my $gene_length = abs ($lend -$rend) + 1;
        foreach my $gene ($self, @isoforms) {
            $gene->{gene_span} = [$lend, $rend];
            $gene->{gene_length} = $gene_length;
        }
    }
    
}




=over 4

=item add_mRNA_exon_obj()

B<Description:> Used to add a single mRNA_exon_obj to the Gene_obj 

B<Parameters:> mRNA_exon_obj

B<Returns:> none

=back

=cut

    ;

sub add_mRNA_exon_obj {
    my ($self) = shift;
    my ($mRNA_exon_obj) = shift;
    if (!ref($self->{mRNA_exon_objs})) {
        $self->{mRNA_exon_objs} = [];
    } 
    my $index = $#{$self->{mRNA_exon_objs}};
    $index++;
    $self->{mRNA_exon_objs}->[$index] = $mRNA_exon_obj;
}

#private
## forcibly set protein sequence value


sub set_protein_sequence {
    my $self = shift;
    my $protein = shift;
    if ($protein) {
        $self->{protein_seq} = $protein;
        $self->{protein_seq_length} = length ($protein);
    } else {
        print STDERR "No incoming protein sequence to set to.\n" . $self->toString();
    }
}

#private
## forcibly set CDS sequence value
sub set_CDS_sequence {
    my $self = shift;
    my $cds_seq = shift;
    if ($cds_seq) {
        $self->{CDS_sequence} = $cds_seq;
        $self->{CDS_sequence_length} = length ($cds_seq);
    } else {
        print STDERR "No incoming CDS sequence to set to\n" . $self->toString();
    }
}

#private
sub set_cDNA_sequence {
    my $self = shift;
    my $cDNA_seq = shift;
    if ($cDNA_seq) {
        $self->{cDNA_sequence} = $cDNA_seq;
        $self->{cDNA_sequence_length} = length($cDNA_seq);
    } else {
        print STDERR "No incoming cDNA sequence to set to.\n" . $self->toString();
    }
}

#private
sub set_gene_sequence {
    my $self = shift;
    my $seq = shift;
    if ($seq) {
        $self->{gene_sequence} = $seq;
        $self->{gene_sequence_length} = length ($seq);
    } else {
        print STDERR "No incoming gene sequence to set to\n" . $self->toString();
    }
}


=over 4

=item create_all_sequence_types()

B<Description:> Given a scalar reference to the genomic sequence, the CDS, cDNA, unspliced transcript and protein sequences are constructed and populated within the Gene_obj

B<Parameters:> $genomic_seq_ref, [%params]

B<Returns:> 0|1

returns 1 upon success, 0 upon failure

By default, the protein and CDS sequence are populated.  If you want the unspliced genomic sequence, you need to specify this in the attributes:

    %params = ( potein => 1,
		CDS => 1,
		cDNA => 1,
		unspliced_transcript => 0)


=back

=cut


## Create all gene sequences (protein, cds, cdna, genomic)
sub create_all_sequence_types {
    my $self = shift;
    my $big_seq_ref = shift;
    my %atts = @_;
    
    unless (ref($big_seq_ref) eq 'SCALAR') {
        print STDERR "I require a sequence reference to create sequence types\n";
        return (undef());
    }
    $self->create_cDNA_sequence($big_seq_ref) unless (exists($atts{cDNA}) && $atts{cDNA});
    $self->create_gene_sequence($big_seq_ref, 1) if ($atts{unspliced_transcript}); #highlight exons by default.
    
    if ($self->is_coding_gene()) {
        $self->create_CDS_sequence ($big_seq_ref) unless (exists ($atts{CDS}) && $atts{CDS});
        $self->create_protein_sequence($big_seq_ref) unless (exists ($atts{protein}) && $atts{protein});
    }
    
    if (my @isoforms = $self->get_additional_isoforms()) {
        foreach my $isoform (@isoforms) {
            $isoform->create_all_sequence_types($big_seq_ref, %atts);
        }
    }
    return(1);
}

#private
## Create cDNA sequence
sub create_cDNA_sequence {
    my $self = shift;
    my $seq_ref = shift;
    my $sequence_ref;
    unless ($seq_ref) {
        print STDERR "The parent sequence must be specified for the cDNA creation method\n";
        return;
    }
    ## hopefully the sequence came in as a reference.  If not, make one to it.
    ## Don't want to pass chromosome sequences in by value!
    if (ref($seq_ref)) {
        $sequence_ref = $seq_ref;
    } else {
        $sequence_ref = \$seq_ref;
    }
    my @exons = $self->get_exons();
    my $strand = $self->{strand};
    my $cDNA_seq = "";
    foreach my $exon_obj (sort {$a->{end5}<=>$b->{end5}} @exons) {
        my $c1 = $exon_obj->{end5};
        my $c2 = $exon_obj->{end3};
        ## sequence retrieval coordinates must be in forward orientation
        my ($coord1, $coord2) = ($strand eq '+') ? ($c1, $c2) : ($c2, $c1);
        $cDNA_seq .= substr ($$sequence_ref, ($coord1 - 1), ($coord2 - $coord1 + 1));
    }
    if ($strand eq '-') {
        $cDNA_seq = &reverse_complement($cDNA_seq);
    }
    $self->set_cDNA_sequence($cDNA_seq);
    return ($cDNA_seq);
}

#private
## create a CDS sequence, and populate the protein field.
sub create_CDS_sequence {
    my $self = shift;
    my $seq_ref = shift;
    my $sequence_ref;
    unless ($seq_ref) {
        print STDERR "The parent sequence must be specified for the CDS creation method\n";
        return;
    }
    
    unless ($self->is_coding_gene()) {
        print STDERR "Warning: No coding region specified for gene: " . $self->toString();
        return("");
    }
    

    ## hopefully the sequence came in as a reference.  If not, make one to it.
    ## Don't want to pass chromosome sequences in by value!
    if (ref($seq_ref)) {
        $sequence_ref = $seq_ref;
    } else {
        $sequence_ref = \$seq_ref;
    }
    my @exons = $self->get_exons();
    my $strand = $self->{strand};
    my $cds_seq = "";
    foreach my $exon_obj (sort {$a->{end5}<=>$b->{end5}} @exons) {
        my $CDS_obj = $exon_obj->get_CDS_obj();
        if (ref $CDS_obj) {
            my ($c1, $c2) = $CDS_obj->get_CDS_end5_end3();
            ## sequence retrieval coordinates must be in forward orientation
            my ($coord1, $coord2) = ($strand eq '+') ? ($c1, $c2) : ($c2, $c1);
            $cds_seq .= substr ($$sequence_ref, ($coord1 - 1), ($coord2 - $coord1 + 1));
        }
    }
    if ($strand eq '-') {
        $cds_seq = &reverse_complement($cds_seq);
    }
    $self->set_CDS_sequence($cds_seq);
    
    return ($cds_seq);
}



sub is_coding_gene {
    my $self = shift;
   
    if ($self->get_CDS_length()) {
        return(1);
    }
    else {
        return(0);
    }
}



#private
## Translation requires parent nucleotide sequence (bac, chromosome, whatever).
sub create_protein_sequence {
    my $self = shift;
    my $seq_ref = shift; # optional
    
    unless ($self->is_coding_gene()) {
        print STDERR "Warning: No coding sequence for gene: " . $self->toString();
        return("");
    }
    
    my $cds_sequence = $self->get_CDS_sequence();
    unless ($cds_sequence) {
            
        ## if has a CDS, then try to translate it if the genome sequence is available.
    
        unless (ref($seq_ref) eq 'SCALAR') {
            print STDERR "I require an assembly sequence ref if the CDS is unavailable\n";
            return;
        }
        $cds_sequence = $self->create_CDS_sequence($seq_ref);
    }
    my $protein = &Nuc_translator::get_protein ($cds_sequence); 
    $self->set_protein_sequence($protein);
    return ($protein);
}

#private
## Create the unspliced nucleotide transcript
sub create_gene_sequence {
    my $self = shift;
    my $big_seq_ref = shift;
    my $highlight_exons_flag = shift; #upcases exons, lowcases introns.
    unless (ref ($big_seq_ref) eq 'SCALAR') {
        print STDERR "I require a reference to the assembly sequence!!\n";
        return (undef());
    }
    my $strand = $self->{strand};
    my ($gene_seq);
    if ($highlight_exons_flag) {
        my @exons = sort {$a->{end5}<=>$b->{end5}} $self->get_exons();
        my $exon = shift @exons;
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        $gene_seq = uc (substr ($$big_seq_ref, $lend - 1, $rend - $lend + 1));
        my $prev_rend = $rend;
        while (@exons) {
            $exon = shift @exons;
            ## Add intron, then exon
            my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
            $gene_seq .= lc (substr ($$big_seq_ref, $prev_rend, $lend - $prev_rend-1));
            $gene_seq .= uc (substr ($$big_seq_ref, $lend - 1, $rend - $lend + 1));
            $prev_rend = $rend;
        }
        
    } else { #just get the sequence spanned by min and max coords
        my ($coord1, $coord2) = sort {$a<=>$b} $self->get_gene_span();
        $gene_seq = substr ($$big_seq_ref, ($coord1 - 1), ($coord2 - $coord1 + 1));
    }
    
    $gene_seq = &reverse_complement($gene_seq) if ($strand eq '-');
    $self->set_gene_sequence($gene_seq);
    return ($gene_seq);
}

## retrieving the sequences

=over 4

=item get_protein_sequence()

B<Description:> Retrieves the protein sequence

B<Parameters:> none

B<Returns:> $protein

Note: You must have called create_all_sequence_types($genomic_ref) before protein sequence is available for retrieval.


=back

=cut
    
    ;

sub get_protein_sequence {
    my $self = shift;
    return ($self->{protein_seq});
}

## alias
sub get_protein_seq {
    my $self = shift;
    return ($self->get_protein_sequence());
}



=over 4

=item get_CDS_sequence()

B<Description:> Retrieves the CDS sequence.  The CDS sequence is the protein-coding nucleotide sequence.

B<Parameters:> none

B<Returns:> $cds

Note: You must have called create_all_sequence_types($genomic_ref) before protein sequence is available for retrieval.

=back

=cut


sub get_CDS_sequence {
    my $self = shift;
    return ($self->{CDS_sequence});
}

=over 4

=item get_cDNA_sequence()

B<Description:> Retrieves the tentative cDNA sequence for the Gene.  The cDNA includes the CDS with potential UTR extensions.

B<Parameters:> none

B<Returns:> $cdna

Note: You must have called create_all_sequence_types($genomic_ref) before protein sequence is available for retrieval.


=back

=cut



sub get_cDNA_sequence {
    my $self = shift;
    return ($self->{cDNA_sequence});
}



sub get_CDS_length {
    my $self = shift;
    
    my $cds_length = 0;
    
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        if (my $cds = $exon->get_CDS_obj()) {
            $cds_length += $cds->length();
        }
    }
    
    
    return ($cds_length);
}

sub get_cDNA_length {
    my $self = shift;
    
    my $cdna_length = 0;
    
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        $cdna_length += $exon->length();
    }

    return($cdna_length);
    
}







=over 4

=item get_gene_sequence()

B<Description:> Retrieves the unspliced transcript of the gene.

B<Parameters:> none

B<Returns:> $unspliced_transcript

=back

=cut


sub get_gene_sequence {
    my $self = shift;
    return ($self->{gene_sequence});
}




=over 4

=item get_gene_synonyms()

B<Description:> Retrieves the Model_feat_name(s) for the synonomous gene models found on other BACs or contigs.

B<Parameters:> none

B<Returns:> @model_feat_names


For Arabidopsis, gene models are found within overlapping regions of BAC sequences, in which the gene models are annotated on both corresponding BACs.  Given a Gene_obj for a model on one BAC, the synomous gene on the overlapping BAC can be identified via this method.


=back

=cut


sub get_gene_synonyms {
    my $self = shift;
    return (@{$self->{gene_synonyms}});
}



=over 4

=item clear_sequence_info()

B<Description:> Clears the sequence fields stored within a Gene_obj, including the CDS, cDNA, gene_sequence, and protein sequence.  Often, these sequence fields, when populated, can consume large amounts of memory in comparison to the coordinate and functional annotation data.  This method is useful to clear this memory when the sequences are not needed.  The create_all_sequence_types($genomic_seq_ref) can be called again later to repopulate these sequences when they are needed.

B<Parameters:> none

B<Returns:> none

=back

=cut


## sequences consume huge amounts of memory in comparison to other gene features.
## want to clear them from time to time to save memory.
    
    ;

sub clear_sequence_info {
    my $self = shift;
    $self->{protein_seq} = undef;   
    $self->{CDS_sequence} = undef;
    $self->{cDNA_sequence} = undef; 
    $self->{gene_sequence} = undef;
}


=over 4

=item set_gene_type()

B<Description:> Sets the type of gene.  Expected types include: 

    protein-coding #default setting
    rRNA
    snoRNA
    snRNA
    tRNA
    
    ...or others as needed.  Nothing is restricted.

B<Parameters:> $type

B<Returns:> none

=back

=cut


####
sub set_gene_type {
    my ($self) = shift;
    my ($gene_type) = shift;
    $self->{gene_type} = $gene_type;
}


=over 4

=item adjust_gene_coordinates()

B<Description:> Used to add or subtract a specified number of bases from each gene component coordinate.

B<Parameters:> $adj_amount

$adj_amoount is a positive or negative integer.

B<Returns:> none

=back

=cut


    ;

####
# add value to all gene component coordinates
sub adjust_gene_coordinates {
    my $self = shift;
    my $adj_amount = shift;
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        my ($end5, $end3) = $exon->get_coords();
        $exon->set_coords($end5 + $adj_amount, $end3 + $adj_amount);
        my $cds = $exon->get_CDS_obj();
        if (ref $cds) {
            my ($end5, $end3) = $cds->get_coords();
            $cds->set_coords($end5 + $adj_amount, $end3 + $adj_amount);
        }
    }
    
    ## don't forget about alt splicing isoforms!
    my @isoforms = $self->get_additional_isoforms();
    foreach my $isoform (@isoforms) {
        $isoform->adjust_gene_coordinates($adj_amount);
    }
    $self->refine_gene_object();
}




=over 4

=item toString()

B<Description:> Textually describes the Gene_obj including coordinates and attributes.

B<Parameters:> <%attributes_list> 

%attributes_list is optional and can control whether certain attributes are included in the textual output

Default settings are:

    %attributes_list = ( 
			 -showIsoforms => 1,  #set to 0 to avoid isoform info to the text output.
			 -showSeqs => 0  #set to 1 for avoiding protein, cds, genomic, cdna seqs as output.
			 )

B<Returns:> $text

=back

=cut

    ;


## retrieve text output describing the gene.
sub toString {
    my $self = shift;
    my %atts = @_;
    # atts defaults:
    #       -showIsoforms=>1
    #       -showSeqs => 0
    
    my $output = "";
    foreach my $key (keys %$self) {
        my $value = $self->{$key};
        unless (defined $value) { next;}
        if (ref $value) {
            if ($key =~ /secondary/ && ref $value eq "ARRAY") {
                foreach my $val (@$value) {
                    $output .= "\t\t$key\t$val\n";
                }
            }
            
            
        } else {
            if ($self->{is_pseudogene} && $key =~ /cds|cdna|protein/i && $key =~ /seq/) {
                next;
            }
            if ((!$atts{-showSeqs}) && $key =~/seq/) { next; }
            if ( ($value eq '0' || !defined($value)) && $key !~/^is_/) { next;} #dont print unpopulated info.
            $output .= "\t$key:\t$value\n";
        }
    }
    $output .= "\tgene_synonyms: @{$self->{gene_synonyms}}\n";
    
    $output .=  "\tmRNA_coords\t";  
    
    if (ref ($self->{mRNA_coords}) eq "HASH") {
        foreach my $end5 (sort {$a<=>$b} keys %{$self->{mRNA_coords}}) {
            $output .=  "$end5-$self->{mRNA_coords}->{$end5} ";
        }
    }
    $output .= "\n"
        . "\tCDS_coords\t";
    if (ref ($self->{CDS_coords}) eq "HASH") {
        foreach my $end5 (sort {$a<=>$b} keys %{$self->{CDS_coords}}) {
            $output .= "$end5-$self->{CDS_coords}->{$end5} ";
        }
    }
    
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        $output .=  "\n\t\tRNA-exon: $exon->{end5}, $exon->{end3}\t";
        my $cds = $exon->{CDS_exon_obj};
        if ($cds) {
            $output .= "CDS-exon: $cds->{end5}, $cds->{end3}";
        }
    }
    
    if (ref $self->{gene_span}) {
        my ($gene_end5, $gene_end3) = @{$self->{gene_span}};
        $output .= "\n\tgene_span: $gene_end5-$gene_end3";
    }
    if (ref $self->{model_span}) {
        my ($model_end5, $model_end3) = @{$self->{model_span}};
        $output .= "\n\tmodel_span: $model_end5-$model_end3";
    }
    my @gene_ontology_objs = $self->get_gene_ontology_objs();
    if (@gene_ontology_objs) {
        $output .= "\n\tGene Ontology Assignments:\n";
        foreach my $go_assignment (@gene_ontology_objs) {
            $output .= "\t" . $go_assignment->toString();
        }
    }
    
    unless (defined ($atts{-showIsoforms}) && $atts{-showIsoforms} == 0) {
        foreach my $isoform ($self->get_additional_isoforms()) {
            $output .= "\n\n\tISOFORM:\n" . $isoform->toString();
        }
    }
    $output .= "\n\n"; #spacer at terminus
    return ($output);
}


####
## Splice site validation section
####

=over 4

=item validate_splice_sites()

B<Description:> Validates the presence of consensus splice sites 

B<Parameters:> $genomic_seq_ref

$genomic_seq_ref is a scalar reference to the string containing the genomic sequence.

B<Returns:> $errors

If the empty string ("") is returned, then no inconsistencies were identified.

=back

=cut

    ;
    
####
sub validate_splice_sites {
    my $self = shift;
    my $asmbl_seq_ref = shift;
    unless (ref ($asmbl_seq_ref)) {
        print STDERR "I require a sequence reference\n";
        return (undef());
    }
    my $error_string = "";
    my $strand = $self->{strand};
    my @exons = $self->get_exons();
    my $num_exons = $#exons + 1;
    if ($num_exons == 1) {
        #no splice sites to confirm.
        return ("");
    }
    for (my $i = 1; $i <= $num_exons; $i++) {
        my $exon_type;
        if ($i == 1) { 
            $exon_type = "initial";
        } elsif ($i == $num_exons) {
            $exon_type = "terminal";
        } else {
            $exon_type = "internal";
        }
        my $exon = $exons[$i - 1]; 
        my ($exon_end5, $exon_end3) = $exon->get_mRNA_exon_end5_end3();
        my ($coord1, $coord2) = sort {$a<=>$b} ($exon_end5, $exon_end3);
        ## get two coordinate sets corresponding to potential splice sites
        my $splice_1_start = $coord1-2-1;
        my $splice_2_start = $coord2-1+1;
        #print "confirming splice sites at "  . ($splice_1_start +1) . " and " . ($splice_2_start + 1) . "\n"if $SEE;
        my $splice_1 = substr ($$asmbl_seq_ref, $splice_1_start, 2);
        my $splice_2 = substr ($$asmbl_seq_ref, $splice_2_start, 2);
        my ($acceptor, $donor) = ($strand eq '+') ? ($splice_1, $splice_2) : (&reverse_complement($splice_2), &reverse_complement($splice_1)); 
        my $check_acceptor = ($acceptor =~ /ag/i);
        my $check_donor = ($donor =~ /gt|gc/i);
        ## associate results of checks with exon type.
        if ($exon_type eq "initial" || $exon_type eq "internal") {
            unless ($check_donor) {
                $error_string .= "non-consensus $donor donor splice site at $coord1\n";
            }
        }
        
        if ($exon_type eq "internal" || $exon_type eq "terminal") {
            unless ($check_acceptor) {
                $error_string .=  "\tnon-consensus $acceptor acceptor splice site at $coord2\n";
            }
        }
    }
    return ($error_string);
}



=over 4

=item get_annot_text()

B<Description:> Provides basic functional annotation for a Gene_obj 

B<Parameters:> none

B<Returns:> $string

$string includes locus, pub_locus, com_name, and pub_comment

=back

=cut


    ;

####
sub get_annot_text {
    my $self = shift;
    my $locus = $self->{locus};
    my $pub_locus = $self->{pub_locus};
    my $com_name = $self->{com_name};
    my $pub_comment = $self->{pub_comment};
    my $text = "";
    foreach my $token ($locus, $pub_locus, $com_name, $pub_comment) {
        if ($token) {
            $text .= "$token ";
        }
    }
    return ($text);
}



=over 4

=item add_isoform()

B<Description:> Adds a Gene_obj to an existing Gene_obj as an alternative splicing variant.

B<Parameters:> Gene_obj

B<Returns:> none

=back

=cut

    ;
sub add_isoform {
    my $self = shift;
    my @gene_objs = @_;
    foreach my $gene_obj (@gene_objs) {
        $self->{num_additional_isoforms}++;
        push (@{$self->{additional_isoforms}}, $gene_obj);
    }
}





=over 4

=item has_additional_isoforms()

B<Description:> Provides number of additional isoforms.  Typically used as a boolean.

B<Parameters:> none

B<Returns:> number of additional isoforms (int)

If no additional isoforms exist, returns 0


boolean usage:

0 = false (has no more)
nonzero = true (has additional isoforms)

=back

=cut

sub has_additional_isoforms {
    my $self = shift;
    return ($self->{num_additional_isoforms});
}



=over 4

=item delete_isoforms()

B<Description:> removes isoforms stored in this Gene_obj (assigning to a new anonymous arrayref)

B<Parameters:> Gene_obj

B<Returns:> none

=back

=cut

sub delete_isoforms {
    my $self = shift;
    $self->{num_additional_isoforms} = 0;
    $self->{additional_isoforms} = [];
}





=over 4

=item get_additional_isoforms()

B<Description:> Retrieves the additional isoforms for a given Gene_obj

B<Parameters:> none

B<Returns:> @Gene_objs

If no additional isoforms exist, an empty array is returned.

=back

=cut


sub get_additional_isoforms {
    my $self = shift;
    return (@{$self->{additional_isoforms}});
}



=over 4

=item get_orientation()

B<Description:> Retrieves the strand orientation of the Gene_obj

B<Parameters:> none

B<Returns:> +|-

=back

=cut


sub get_orientation {
    my $self = shift;
    return ($self->{strand});
}



sub get_strand { ## preferred
	my $self = shift;
	return($self->get_orientation());
}



=over 4

=item add_gene_ontology_objs()

B<Description:> Adds a list of Gene_ontology objects to a Gene_obj

B<Parameters:> @Gene_ontology_objs

@Gene_ontology_objs is a list of objects instantiated from Gene_ontology.pm

B<Returns:> none

=back

=cut


sub add_gene_ontology_objs {
    my ($self, @ontology_objs) = @_;
    push (@{$self->{GeneOntology}}, @ontology_objs);
}



=over 4

=item get_gene_ontology_objs()

B<Description:> Retrieves Gene_ontology objs assigned to the Gene_obj

B<Parameters:> none

B<Returns:> @Gene_ontology_objs

@Gene_ontology_objs are objects instantiated from package Gene_ontology  (See Gene_ontology.pm)

=back

=cut

    ;

sub get_gene_ontology_objs {
    my $self = shift;
    if (ref ($self->{GeneOntology})) {
        return (@{$self->{GeneOntology}});
    } else {
        return (());
    }
}


=over 4

=item set_5prime_partial()

B<Description:> Sets the status of the is_5prime_partial attribute

B<Parameters:> 1|0

B<Returns:> none


5prime partials are partial on their 5prime end and lack start codons.


=back

=cut

sub set_5prime_partial() {
    my $self = shift;
    my $value = shift;
    $self->{is_5prime_partial} = $value;
}



=over 4

=item set_3prime_partial()

B<Description:> Sets the is_3prime_partial status

B<Parameters:> 1|0

B<Returns:> none

3prime partials are partial on their 3prime end and lack stop codons.

=back

=cut


sub set_3prime_partial() {
    my $self = shift;
    my $value = shift;
    $self->{is_3prime_partial} = $value;
}



=over 4

=item is_5prime_partial()

B<Description:> Retrieves the 5-prime partial status of the gene.

B<Parameters:> none

B<Returns:> 1|0

=back

=cut


sub is_5prime_partial() {
    my $self = shift;
    return ($self->{is_5prime_partial});
}


=over 4

=item is_3prime_partial()

B<Description:> Retrieves the 3-prime partial status of the gene.

B<Parameters:> none

B<Returns:> 1|0

=back

=cut


sub is_3prime_partial() {
    my $self = shift;
    return ($self->{is_3prime_partial});
}

=over 4

=item get_5prime_UTR_coords
	

B<Description:> returns a list of coordinate pairs corresponding to the 5\' UTR coordinates

B<Parameters:> none

B<Returns:> ([end5,end3], ...) or empty list if none exist

=back

=cut


    ;

sub get_5prime_UTR_coords {
    my $self = shift;
    
    my $strand = $self->get_orientation();
    
    my @exons = $self->get_exons();
    
    my $seen_CDS_flag = 0;
    
    my @utr_coords;
    foreach my $exon (@exons) { #relying on a sorted list
        my ($exon_end5, $exon_end3) = $exon->get_coords();
        if (my $cds = $exon->get_CDS_obj()) {
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            if ($exon_end5 != $cds_end5) {
                my $adj_utr_end3_coord = ($strand eq '+') ? ($cds_end5 -1) : ($cds_end5 +1);
                push (@utr_coords, [$exon_end5, $adj_utr_end3_coord]);
            } 
            
            $seen_CDS_flag = 1;
            
        } else {
            push (@utr_coords, [$exon_end5, $exon_end3]);
        }
        
        if ($seen_CDS_flag) {
            last;
        }
        
    }
    
    return (@utr_coords);
}



=over 4

=item get_3prime_UTR_coords
	

B<Description:> returns a list of coordinate pairs corresponding to the 3\' UTR coordinates

B<Parameters:> none

B<Returns:> ([end5,end3], ...) or empty list if none exist

=back

=cut

    ;

sub get_3prime_UTR_coords {
    my $self = shift;
    
    my $strand = $self->get_orientation();
    
    my @exons = reverse $self->get_exons();
    
    my @utr_coords;
    my $seen_CDS_flag = 0;
    foreach my $exon (@exons) { #relying on a reverse sorted list (3' exons should come first)
        my ($exon_end5, $exon_end3) = $exon->get_coords();
        if (my $cds = $exon->get_CDS_obj()) {
            $seen_CDS_flag = 1;
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            if ($exon_end3 != $cds_end3) {
                my $adj_utr_end5_coord = ($strand eq '+') ? ($cds_end3 +1) : ($cds_end3 -1);
                push (@utr_coords, [$adj_utr_end5_coord, $exon_end3]);
            } 
	  	    
        } else {
            push (@utr_coords, [$exon_end5, $exon_end3]);
        }
        if ($seen_CDS_flag) { 
            last;
        }
    }

    if (@utr_coords) {
        @utr_coords = reverse @utr_coords;
    }
    
    return (@utr_coords);
}





=over 4

=item has_UTRs()

B<Description:> indicates presence of UTR annotated in Gene 

B<Parameters:> none

B<Returns:> ( has_5prime_UTR() || has_3prime_UTR() )

=back

=cut

sub has_UTRs {
    my $self = shift;
    return ( ($self->has_5prime_UTR() || $self->has_3prime_UTR() ) );
}



####
sub has_5prime_UTR {
    my $self = shift;
    return ($self->get_5prime_UTR_length() > 2);
}

####
sub has_3prime_UTR {
    my $self = shift;
    return($self->get_3prime_UTR_length() > 2);
}


####
sub get_5prime_UTR_length {
    my $self = shift;
    my @prime5_UTR_coords = $self->get_5prime_UTR_coords();
    
    my $len = 0;
    for my $coordset (@prime5_UTR_coords) {
        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
        $len += ($rend - $lend) + 1;
    }
    return($len);
}


####
sub get_3prime_UTR_length {
    my $self = shift;
    my @prime3_UTR_coords = $self->get_3prime_UTR_coords();
    
    my $len = 0;
    for my $coordset (@prime3_UTR_coords) {
        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
        $len += ($rend - $lend) + 1;
    }
    return($len);
}
    


=over 4

=item get_5prime_UTR_sequence()

B<Description:> retrieves 5prime UTR sequence

B<Parameters:> genome sequence reference

B<Returns:> string

=back

=cut

####
sub get_5prime_UTR_sequence {
    my $self = shift;
    my ($genome_seq_ref) = @_;
    unless (ref $genome_seq_ref eq "SCALAR") {
        confess "error, require genome sequence string reference";
    }

    unless ($self->has_5prime_UTR()) {
        return "";
    }
    
    my $orientation = $self->get_orientation();
    my @coords = $self->get_5prime_UTR_coords();
    
    @coords = sort {$a->[0]<=>$b->[0]} @coords;

    my $UTR_seq = "";
    foreach my $coordset (@coords) {
        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
        
        my $length = $rend - $lend + 1;
        $UTR_seq .= substr($$genome_seq_ref, $lend - 1, $length);
    }

    if ($orientation eq '-') {
        $UTR_seq = &reverse_complement($UTR_seq);
    }

    ## verify:
    $self->create_all_sequence_types($genome_seq_ref);
    my $cDNA = $self->get_cDNA_sequence();
    

    unless (index($cDNA, $UTR_seq) == 0) {
        confess "Error, couldn't find UTR in cDNA";
    }

    
    return ($UTR_seq);
}
        

=over 4

=item get_3prime_UTR_sequence()

B<Description:> retrieves 5prime UTR sequence

B<Parameters:> genome sequence reference

B<Returns:> string

=back

=cut

####
sub get_3prime_UTR_sequence {
    my $self = shift;
    my ($genome_seq_ref) = @_;
    unless (ref $genome_seq_ref eq "SCALAR") {
        confess "error, require genome sequence string reference";
    }

    unless ($self->has_3prime_UTR()) {
        return "";
    }
    
    my $orientation = $self->get_orientation();
    my @coords = $self->get_3prime_UTR_coords();
    
    @coords = sort {$a->[0]<=>$b->[0]} @coords;

    my $UTR_seq = "";
    foreach my $coordset (@coords) {
        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
        
        my $length = $rend - $lend + 1;
        $UTR_seq .= substr($$genome_seq_ref, $lend - 1, $length);
    }

    if ($orientation eq '-') {
        $UTR_seq = &reverse_complement($UTR_seq);
    }

    ## verify:
    $self->create_all_sequence_types($genome_seq_ref);
    my $cDNA = $self->get_cDNA_sequence();
    my $cDNA_length = length($cDNA);
    my $utr_length = length($UTR_seq);

    my $utr_start_pos = $cDNA_length - $utr_length + 1;

    unless ((my $cDNA_utr =  lc substr($cDNA, $utr_start_pos - 1, $utr_length)) eq lc $UTR_seq) {
        confess "Error, 3' UTR extracted from cDNA is different from UTR sequence extracted from genome.\n"
            . "cDNA_utr:\n$cDNA_utr\nUTR_from_genome:\n$UTR_seq\n\n";
    }
    
    
    return ($UTR_seq);
}
        



=over 4

=item trim_UTRs()

B<Description:> Trims the UTR of the Gene_obj so that the Exon coordinates are identical to the CDS coordinates.  Exons which lack CDS components and are completely UTR are removed. 

B<Parameters:> none

B<Returns:> none

=back

=cut

    ;

sub trim_UTRs {
    my $self = shift;
    
    ## adjust exon coordinates to CDS coordinates.
    ## if cds doesn't exist, rid exon:
    
    my @new_exons;
    
    my @exons = $self->get_exons();
    foreach my $exon (@exons) {
        if (my $cds = $exon->get_CDS_obj()) {
            my ($exon_end5, $exon_end3) = $exon->get_coords();
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            
            if ($exon_end5 != $cds_end5 || $exon_end3 != $cds_end3) {
                $exon->set_coords($cds_end5, $cds_end3);
            }
            push (@new_exons, $exon);
        }
    }
    $self->{mRNA_exon_objs} = 0; #clear current gene structure
    $self->{mRNA_exon_objs} = \@new_exons; #replace gene structure
    $self->refine_gene_object(); #update
    return ($self);
}




=over 4
    
=item remove_CDS_exon()

B<Description:> Removes any existing CDS_exon_obj from this mRNA_exon_obj

B<Parameters:> none

B<Returns:> none

=back

=cut

sub remove_CDS_exon {
    my $self = shift;
    $self->{CDS_exon_obj} = 0;
}





=over 4

=item get_gene_names()

B<Description:> Retrieves gene names  (primary gene name followed by secondary gene names, "$;" delimited.

B<Parameters:> none

B<Returns:> string
				       
     see $gene_obj->{gene_name}
     see $gene_obj->get_secondary_names()

secondary gene names sorted lexicographically


=back

=cut




####
sub get_gene_names {
    my $gene_obj = shift;
    my @gene_names;
    if ($gene_obj->{gene_name}) {
        push (@gene_names, $gene_obj->{gene_name});
    }
    if (my @secondary_names = $gene_obj->get_secondary_gene_names()) {
        push (@gene_names, @secondary_names);
    }
    my $ret_gene_names = join ("$;" , @gene_names);
    return ($ret_gene_names);
}



=over 4

=item get_secondary_gene_names()

B<Description:> Retrieves secondary gene names as a "$;" delimited string.

B<Parameters:> none

B<Returns:> string

=back

=cut


####
sub get_secondary_gene_names {
    my ($gene_obj) = @_;
    return (sort @{$gene_obj->{secondary_gene_names}});
}




=over 4

=item get_product_names()

B<Description:> Retrieves product name, with the primary product name followed by secondary product names, delimited by "$;"

B<Parameters:> none

B<Returns:> string

    see $gene_obj->{com_name} for primary product name
    see $gene_obj->get_secondary_product_names()

=back

=cut

    ;

####
sub get_product_names {
    my $gene_obj = shift;
    my @product_names;
    if ($gene_obj->{com_name}) {
        push (@product_names, $gene_obj->{com_name});
    }
    if (my @secondary_names = $gene_obj->get_secondary_product_names()) {
        push (@product_names, @secondary_names);
    }
    my $ret_product_names = join ("$;", @product_names);
    return ($ret_product_names);
}



=over 4

=item get_secondary_product_names()

B<Description:> Retrieves secondary product names, delimited by "$;" and sorted lexicographically.

B<Parameters:> none 

B<Returns:> string

=back

=cut


####
sub get_secondary_product_names {
    my ($gene_obj) = @_;
    return (sort @{$gene_obj->{secondary_product_names}});
}



=over 4

=item get_gene_symbols()

B<Description:> Retrieves primary gene symbol followed by secondary gene symbols, delimited by "$;"

B<Parameters:> none

B<Returns:> string

    see $gene_obj->{gene_sym}
    see $gene_obj->get_secondary_gene_symbols()

=back

=cut

    ;

####
sub get_gene_symbols {
    my $gene_obj = shift;
    my @gene_symbols;
    if ($gene_obj->{gene_sym}) {
        push (@gene_symbols, $gene_obj->{gene_sym});
    }
    if (my @secondary_symbols = $gene_obj->get_secondary_gene_symbols()) {
        push (@gene_symbols, @secondary_symbols);
    }
    my $ret_gene_symbols = join ("$;", @gene_symbols);
    return ($ret_gene_symbols);
}


=over 4

=item get_secondary_gene_symbols()

B<Description:> Retrieves secondary gene symbols, delimited by "$;" and sorted lexicographically

B<Parameters:> none

B<Returns:> string

=back

=cut


####
sub get_secondary_gene_symbols {
    my ($gene_obj) = @_;
    return (sort @{$gene_obj->{secondary_gene_symbols}});
}



=over 4

=item get_ec_numbers()

B<Description:> Retrieves primary EC number followed by secondary EC numbers, "$;" delimited

B<Parameters:> none

B<Returns:> string

    see $gene_obj->{ec_num}
    see $gene_obj->get_secondary_ec_numbers()
    
=back

=cut

    ;

####
sub get_ec_numbers {
    my $gene_obj = shift;
    my @ec_numbers;
    if ($gene_obj->{ec_num}) {
        push (@ec_numbers, $gene_obj->{ec_num});
    }
    if (my @secondary_ec_numbers = $gene_obj->get_secondary_ec_numbers()) {
        push (@ec_numbers, @secondary_ec_numbers);
    }
    my $ret_ec_numbers = join ("$;", @ec_numbers);
    return ($ret_ec_numbers);
}



=over 4

=item get_secondary_ec_numbers()

B<Description:> Retrieves secondary EC numbers, "$;" delimited and sorted lexicographically

B<Parameters:> none

B<Returns:> string


=back

=cut


####
sub get_secondary_ec_numbers {
    my ($gene_obj) = @_;
    return (sort @{$gene_obj->{secondary_ec_numbers}});
}



=over 4

=item add_secondary_gene_names()

B<Description:> Adds secondary gene name(s) 

B<Parameters:> (gene_name_1, gene_name_2, ....)

Single gene name or list of gene names is allowed


B<Returns:> none

=back

=cut



####
sub add_secondary_gene_names {
    my ($gene_obj, @gene_names) = @_;
    push (@{$gene_obj->{secondary_gene_names}}, @gene_names);
}


=over 4

=item add_secondary_product_names()

B<Description:> Adds secondary product names

B<Parameters:> (product_name_1, product_name_2, ...)

Single or list of product names as parameter

B<Returns:> none

Primary gene name added directly as an attribute like so
    $gene_obj->{gene_name} = name

=back

=cut


####
sub add_secondary_product_names {
    my ($gene_obj, @product_names) = @_;
    &trim_leading_trailing_ws(\@product_names);
    push (@{$gene_obj->{secondary_product_names}}, @product_names);
}


=over 4

=item add_secondary_gene_symbols()

B<Description:> Add secondary gene symbols

B<Parameters:> (gene_symbol_1, gene_symbol_2, ...)

String or list context

B<Returns:> none

Primary gene_symbol added directly as attribute like so:
    $gene_obj->{gene_sym} = symbol

=back

=cut


####
sub add_secondary_gene_symbols {
    my ($gene_obj, @gene_symbols) = @_;
    &trim_leading_trailing_ws(\@gene_symbols);
    push (@{$gene_obj->{secondary_gene_symbols}}, @gene_symbols);
}





=over 4

=item add_secondary_ec_numbers()

B<Description:> Add secondary Enzyme Commission (EC) numbers

B<Parameters:> (EC_1, EC_2, ...)

String or list context

B<Returns:> none


Primary EC number added directly as an attribute like so:
    $gene_obj->{ec_num} = EC_number

=back

=cut


####
sub add_secondary_ec_numbers {
    my ($gene_obj, @ec_numbers) = @_;
    &trim_leading_trailing_ws(\@ec_numbers);
    push (@{$gene_obj->{secondary_ec_numbers}}, @ec_numbers);
}

####
sub to_alignment_GFF3_format {
    my ($gene_obj, $id, $target, $source) = @_;

    unless (defined $source) {
        $source = ".";
    }

    ## Note, only examines gene_obj and doesn't go deeper into alt-splicing layers, ... send isoforms in as separate objs.

    unless ( (ref $gene_obj)  && defined($id) && defined($target)) {
        croak "Error, need gene_obj, id, and target names as params";
    }

    my $gff3_alignment_text = "";
        
    my $orient = $gene_obj->get_orientation();
    my $scaff = $gene_obj->{asmbl_id};
    
    my @exons = sort {$a->{end5}<=>$b->{end5}} $gene_obj->get_exons();
    
    if ($orient eq '-') {
        @exons = reverse @exons;
    }
    

    my $match_lend = 0;

    foreach my $exon (@exons) {
        
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        
        my $m_lend = $match_lend + 1;
        my $m_rend = $match_lend + ($rend - $lend + 1);
        

        $gff3_alignment_text .= join("\t", $scaff, $source, "match", $lend, $rend, "100", $orient, '.', # giving everything 100% identity since genome-based 
                                     "ID=$id;Target=$target $m_lend $m_rend +") . "\n";
        
        
        $match_lend = $m_rend;
        
        
    }
    
    return($gff3_alignment_text);
}




####
sub to_transcript_GTF_format {
	my ($gene_obj) = @_;

	## no worries about protein-coding regions.  Only report transcripts and exons tied to a particular gene.
	## used with cufflinks package for computing FPKM values

	my $gtf_text = "";

	foreach my $gene ($gene_obj, $gene_obj->get_additional_isoforms()) {
		
		my $gene_id = $gene->{TU_feat_name} || "";
		my $transcript_id = $gene->{Model_feat_name} || "";
		my $asmbl_id = $gene_obj->{asmbl_id};
		my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_transcript_span();
		my $orientation = $gene_obj->get_orientation();
        
        my $com_name = $gene_obj->{com_name} || "";
        $com_name =~ s/;/_/g;
        $com_name =~ s/\"//g;
        
        
        if ($gene->{gene_type} eq "protein-coding") {
            my @exons = $gene->get_exons();
            
            $gtf_text .= join("\t", $asmbl_id, ".", "transcript", $lend, $rend, ".", $orientation, ".", 
                              "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; name \"$com_name\";") . "\n";
            
            foreach my $exon (@exons) {
                my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
                
                $gtf_text .= join("\t", $asmbl_id, ".", "exon", $lend, $rend, ".", $orientation, ".", 
                                  "gene_id \"$gene_id\"; transcript_id \"$transcript_id\";") . "\n";
                
            }
            
        }
        else {
            
            ## non-protein-coding features
            $gtf_text .= join("\t", $asmbl_id, ".", $gene->{gene_type}, $lend, $rend, ".", $orientation, ".", 
                                  "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; name \"$com_name\";") . "\n";
            
            
        }
        
        $gtf_text .= "\n";
    }
	
	
	return($gtf_text);
}



=over 4

=item to_GTF_format()

B<Description:> Outputs text corresponding to the representation of the gene in GTF format.

B<Parameters:> $genome_seq_ref,  %preferences

B<Returns:> string


GTF format is described in "Current Protocols in Bioinformatics(2003)" 4.8.1-4.8.19
in "Using TWINSCAN to Predict Gene Structures in Genomic DNA Sequences".

Each line of the GTF format includes the following tab-delimited fields:

[seqname] [source] [feature] [start] [end] [score] [strand] [frame] [attributes]

This is further elaborated below:

[feature] contains one of the following: start_codon, stop_codon, CDS
[attributes] contains 'gene_id' and 'transcript_id' fields.  All features of the same transcript should share the same transcript_id value.  By default, the TU_feat_name and model_feat_name are used as the gene_id and transcript_id, respectively.


    Using the %preferences input parameter, the preferred values or gene attributes can be used for seqname, source, gene_id, or transcript_id, each used as a key to the %preferences hash.  Given the value of %preferences is a gene attribute, that attribute value will be used, otherwise, the raw value will be used.

For example:  %preferences = ( seqname => 'mySeqname',
                               gene_id => 'pub_locus' );

Would result in 'mySeqname' used in the [seqname] field, and the $gene_obj->{pub_locus} value  

Here are the defaults:
[seqname] = asmbl_id
[source] = annotation
gene_id (TU_feat_name)
transcript_id (Model_feat_name)

** Partial Genes are NOT Supported **  ( undef is returned )
** Genes with split start or stop codons are unsupported ** (undef is returned)

=back

=cut

    ;

sub to_GTF_format {
    my $gene_obj = shift;
    my ($genome_seq_ref, %preferences) = @_;
    
    unless (ref $genome_seq_ref) { 
        confess "Error, need genome seq reference as param";
    }
    
    my $is_pseudogene = $gene_obj->is_pseudogene();
    

    my $TU_feat_name = $gene_obj->{TU_feat_name};
    my $model_feat_name = $gene_obj->{Model_feat_name};
    
    # rid whitespace in identifiers
    $TU_feat_name =~ s/\s+/_/g;
    $model_feat_name =~ s/\s+/_/g;
    
    my $seqname = $preferences{seqname} || $gene_obj->{asmbl_id};
    my $source = $preferences{source} || $gene_obj->{source} || ".";
    
    my $gene_id;
    if (my $token = $preferences{gene_id}) {
        $gene_id = $gene_obj->{$token};
    } else {
        $gene_id = $TU_feat_name;
    }
    
    my $transcript_id;
    if (my $token = $preferences{model_id}) {
        $transcript_id = $gene_obj->{$token};
    } else {
        $transcript_id = $model_feat_name;
    }
        
    my @exons = $gene_obj->get_exons();
    my $orientation = $gene_obj->get_orientation();    
    my @gtf_text;
    
    my $gene_obj_for_gtf = $gene_obj;  #if got stop codon, will need to strip it off.
	my $com_name = $gene_obj->{com_name};
	$com_name =~ s/\s+$// if $com_name;
	$com_name =~ s/[\"\']//g if $com_name;
	
	my $name_txt = ($com_name) ? "Name \"$com_name\";" : "";


    ## Gene record
    unless ($preferences{'gene_record_already_done'}) {
        
        my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj->get_gene_span();
        
        push (@gtf_text, [$seqname,
                          $source,
                          "gene",
                          $gene_lend,
                          $gene_rend,
                          "0",
                          $orientation,
                          ".",
                          "gene_id \"$gene_id\"; $name_txt"]);
    }
    
    ## Transcript record
    my ($trans_lend, $trans_rend) = sort {$a<=>$b} $gene_obj->get_transcript_span();
    push (@gtf_text, [$seqname,
                      $source,
                      "transcript",
                      $trans_lend,
                      $trans_rend,
                      "0",
                      $orientation,
                      ".",
                      "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; $name_txt"]);
    
    unless ($is_pseudogene) {
        $gene_obj->set_CDS_phases($genome_seq_ref);
    
        
        ## check for start and stop codons.
        my $cds_seq = uc $gene_obj->create_CDS_sequence($genome_seq_ref);
        my @stop_codons = &Nuc_translator::get_stop_codons();
        
        my $first_CDS_segment = $gene_obj->get_first_CDS_segment();
        my $first_phase = $first_CDS_segment->get_phase();
        my $cds_is_integral_codon_num = (length($cds_seq) % 3 == 0) ? 1 : 0;
        
        ## examine start codon:
        my $init_codon = substr($cds_seq, 0, 3);
        if ($first_phase == 0 && $init_codon eq 'ATG') { # got start codon.
            my @start_coordsets = $gene_obj->get_start_codon_coordinates();
            foreach my $start_pair (@start_coordsets) {
                my ($start_lend, $start_rend) = sort {$a<=>$b} @$start_pair;
                push (@gtf_text, [$seqname,
                                  $source,
                                  "start_codon",
                                  $start_lend,
                                  $start_rend,
                                  "0",
                                  $orientation,
                                  "0",
                                  "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; $name_txt"]);
            }
        }
        
        my $candidate_stop_codon = uc substr($cds_seq, length($cds_seq) - 3, 3);
        my @found_stop = grep { $_ eq $candidate_stop_codon } @stop_codons;
        
        if (@found_stop) {
            # got a stop codon.
            # check to see that the stop codon is in-frame.
            if ((length($cds_seq) - $first_phase) % 3 == 0) { # yes, stop is in frame.
                
                my @stop_codon_coords = $gene_obj->get_stop_codon_coords();
                foreach my $stop_pair (@stop_codon_coords) {
                    my ($stop_lend, $stop_rend) = sort {$a<=>$b} @$stop_pair;
                    
                    push (@gtf_text, [$seqname,
                                      $source,
                                      "stop_codon",
                                      $stop_lend,
                                      $stop_rend,
                                      "0",
                                      $orientation,
                                      "0",
                                      "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; $name_txt"]);
                }
                
                $gene_obj_for_gtf = $gene_obj->clone_gene();
                
                $gene_obj_for_gtf->trim_stop_codon();
                
            }
        }
    
    }
    
    ## report the exons and CDS regions:
    foreach my $exon ($gene_obj_for_gtf->get_exons()) {

        my ($exon_lend, $exon_rend) = sort {$a<=>$b} $exon->get_coords();
        
        push (@gtf_text, [$seqname,
                          $source,
                          "exon",
                          $exon_lend,
                          $exon_rend,
                          "0",
                          $orientation,
                          ".",
                          "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; $name_txt"]);
        

        
        my $cds = ($is_pseudogene) ? $exon : $exon->get_CDS_exon_obj();
        
        if ($cds) {
            my $phase = ".";
            unless ($is_pseudogene) {
                $phase = $cds->get_phase();
                if ($phase) {
                    $phase = ($phase == 1) ? 2 : 1; # reverse it according to GFF3 vs. GTF representation.
                }
            }
            
            my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds->get_coords();
         
            push (@gtf_text, [$seqname,
                              $source,
                              "CDS",
                              $cds_lend,
                              $cds_rend,
                              "0",
                              $orientation,
                              "$phase",
                              "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; $name_txt"]);
        }
        
            
    }
    
    unless ($is_pseudogene) {

        ## Get UTR info:
        {
            for my $pair ($gene_obj->get_3prime_UTR_coords) {
                my ($lend,$rend) = sort {$a<=>$b} @$pair;
                push (@gtf_text,   [$seqname,
                                    $source,
                                    "3UTR",
                                    $lend,
                                    $rend,
                                    "0",
                                    $orientation,
                                    "0",
                                    "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; $name_txt"] );
            }
            for my $pair ($gene_obj->get_5prime_UTR_coords) {
                my ($lend,$rend) = sort {$a<=>$b} @$pair;
                push (@gtf_text,   [$seqname,
                                    $source,
                                    "5UTR",
                                    $lend, 
                                    $rend,
                                    "0",
                                    $orientation,
                                    "0",
                                    "gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; $name_txt" ] );
            }
        }
        
    }
    
    @gtf_text = sort {$a->[3] <=> $b->[3]} @gtf_text;
    
    if ($orientation eq '-') {
        @gtf_text = reverse @gtf_text;
    }
        
    my $GTF = "";
    foreach my $gtf_row (@gtf_text) {
        $GTF .= join ("\t", @$gtf_row) . "\n";
    }
    
    foreach my $isoform ($gene_obj->get_additional_isoforms()) {
        my %iso_pref = %preferences;
        $iso_pref{'gene_record_already_done'} = 1;
        $GTF .= "\n" . $isoform->to_GTF_format($genome_seq_ref, %iso_pref);
    }
    
    return ($GTF);
}





####
sub get_start_codon_coordinates {
    my $gene_obj = shift;
    
    my $orient = $gene_obj->get_orientation();
    
    ## just want the coordinate pairs that define the first three CDS bases.
    
    my @cds_coords;
    foreach my $exon ($gene_obj->get_exons()) {
        if (my $cds = $exon->get_CDS_exon_obj()) {
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            push (@cds_coords, [$cds_end5, $cds_end3]);
        }
    }
    
    my @start_coords;
    my $start_len_want = 3;
    foreach my $cds_coordpair (@cds_coords) {
        my ($cds_end5, $cds_end3) = @$cds_coordpair;
        my $cds_seg_len = abs ($cds_end3 - $cds_end5) + 1;
        
        my $extract_len = ($cds_seg_len < $start_len_want) ? $cds_seg_len : $start_len_want;
        if ($orient eq '+') {
            push (@start_coords, [$cds_end5, $cds_end5 + $extract_len - 1]);
        }
        else {
            push (@start_coords, [$cds_end5, $cds_end5 - $extract_len + 1]);
        }
        $start_len_want -= $extract_len;
        
        if ($start_len_want <= 0) { last; }
    }
    
    if ($start_len_want > 0) { 
        confess "Error, trouble extracting start codon coordinates from cds coordsets: " . Dumper (\@cds_coords);
    }
    
    return (@start_coords);
}




####
sub get_stop_codon_coords {
    my $gene_obj = shift;
    
    my $orient = $gene_obj->get_orientation();
    
    ## just want the coordinate pairs that define the last three CDS bases.
    
    my @cds_coords;
    foreach my $exon (reverse $gene_obj->get_exons()) {
        if (my $cds = $exon->get_CDS_exon_obj()) {
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            push (@cds_coords, [$cds_end5, $cds_end3]);
        }
    }
    
    my @stop_coords;
    my $stop_len_want = 3;
    foreach my $cds_coordpair (@cds_coords) {
        my ($cds_end5, $cds_end3) = @$cds_coordpair;
        my $cds_seg_len = abs ($cds_end3 - $cds_end5) + 1;
        
        my $extract_len = ($cds_seg_len < $stop_len_want) ? $cds_seg_len : $stop_len_want;
        if ($orient eq '+') {
            push (@stop_coords, [$cds_end3 - $extract_len + 1, $cds_end3]);
        }
        else {
            push (@stop_coords, [$cds_end3, $cds_end3 + $extract_len - 1]);
        }
        $stop_len_want -= $extract_len;
        
        if ($stop_len_want <= 0) { last; }
    }
    
    if ($stop_len_want > 0) { 
        confess "Error, trouble extracting stop codon coordinates from cds coordsets: " . Dumper (\@cds_coords);
    }
    
    
    return (@stop_coords);
    
}


####
sub trim_stop_codon {
    my $gene_obj = shift;
    
    ## just trimming the last three bases from the CDS's, changing the current gene object.
    
    my @exons = reverse $gene_obj->get_exons();
    
    my $orient = $gene_obj->get_orientation();
    
    my $stop_len_want = 3;
    foreach my $exon (@exons) {
        if (my $cds = $exon->get_CDS_exon_obj()) {
            
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            my $cds_seg_len = abs ($cds_end3 - $cds_end5) + 1;
            
            my $extract_len = ($cds_seg_len < $stop_len_want) ? $cds_seg_len : $stop_len_want;
            
            if ($cds_seg_len == $extract_len) {
                # delete it!
                $exon->delete_CDS_exon_obj();
            }
            else {
                ## truncate it by extract_len
                if ($orient eq '+') {
                    $cds->{end3} -= $extract_len;
                }
                
                else {
                    $cds->{end3} += $extract_len;
                }
            }
            $stop_len_want -= $extract_len;
            
            if ($stop_len_want <= 0) { last; }
        }
    }
    if ($stop_len_want > 0) { 
        confess "Error, trouble extracting all stop codon coordinates from cds coordsets. " . $gene_obj->toString();
    }
    
    return;
    
}



=over 4
    
=item to_GFF3_format()

B<Description:> Outputs text corresponding to the representation of the gene in GFF3 format (still under development).

B<Parameters:> 

B<Returns:> string

GFF3 defined at:
http://song.sourceforge.net/gff3-jan04.shtml

(some text lifted from above site provided below for reference purposes)

The format consists of 9 columns, separated by tabs or spaces.  The
following unescaped characters are allowed within fields:
[a-zA-Z0-9.:^*$@!+_?-].  All other characters must must be escaped
using the URL escaping conventions.  Unescaped quotation marks,
backslashes and other ad-hoc escaping conventions that have been added
to the GFF format are explicitly forbidden.  The =, ; and % characters
have reserved meanings as described below, and must be escaped when
used in other contexts.

Undefined fields are replaced with the "." character, as described in
the original GFF spec.

Column 1: "seqid"

The ID of the landmark used to establish the coordinate system for the
current feature.  IDs must contain alphanumeric characters.
Whitespace, if present, must be escaped using URL escaping rule
(e.g. space="%20" or "+").  Sequences must *NOT* begin with an
unescaped ">".

Column 2: "source"

The source of the feature.  This is unchanged from the older GFF specs
and is not part of a controlled vocabulary.

Column 3: "type"

The type of the feature (previously called the "method").  This is
constrained to be either: (a) a term from the "lite" sequence
ontology, SOFA; or (b) a SOFA accession number.  The latter
alternative is distinguished using the syntax SO:000000.

Columns 4 & 5: "start" and "end"

The start and end of the feature, in 1-based integer coordinates,
relative to the landmark given in column 1.  Start is always less than
or equal to end.

For zero-length features, such as insertion sites, start equals end
and the implied site is to the right of the indicated base.  This
convention holds regardless of the strandedness of the feature.

Column 6: "score"

The score of the feature, a floating point number.  As in earlier
versions of the format, the semantics of the score are ill-defined.
It is strongly recommended that E-values be used for sequence
similarity features, and that P-values be used for ab initio gene
prediction features.

Column 7: "strand"

The strand of the feature.  + for positive strand (relative to the
landmark), - for minus strand, and . for features that are not
stranded.  In addition, ? can be used for features whose strandedness
is relevant, but unknown.

Column 8: "phase"

For features of type "exon", the phase indicates where the feature
begins with reference to the reading frame.  The phase is one of the
integers 0, 1,or 2, indicating that the first base of the feature
corresponds to the first, second or last base of the codon,
respectively.  This is NOT to be confused with the frame, but relates
to the relative position of the translational start in whatever strand
the feature is in.

Column 9: "attributes"

A list of feature attributes in the format tag=value.  Multiple
tag=value pairs are separated by semicolons.  URL escaping rules are
used for tags or values containing the following characters: ",=;".
Whitespace should be replaced with the "+" character or the %20 URL
escape.  This will allow the file to survive text processing programs
that convert tabs into spaces.

These tags have predefined meanings:

    ID	   Indicates the name of the feature.  IDs must be unique
	   within the scope of the GFF file.

    Name   Display name for the feature.  This is the name to be
           displayed to the user.  Unlike IDs, there is no requirement
	   that the Name be unique within the file.

    Alias  A secondary name for the feature.  It is suggested that
	   this tag be used whenever a secondary identifier for the
	   feature is needed, such as locus names and
	   accession numbers.  Unlike ID, there is no requirement
	   that Alias be unique within the file.

    Parent Indicates the parent of the feature.  A parent ID can be
	   used to group exons into transcripts, transcripts into
	   genes, an so forth.  A feature may have multiple parents.

    Target Indicates the target of a nucleotide-to-nucleotide or
	   protein-to-nucleotide alignment.  The format of the
	   value is "target_id+start+end".

    Gap    The alignment of the feature to the target if the two are
          not colinear (e.g. contain gaps).  The alignment format is
	  taken from the CIGAR format described in the
	  Exonerate documentation.
	  (http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate
           ?cvsroot=Ensembl).  See "THE GAP ATTRIBUTE" for a description
	   of this format.

    Note   A free text note.

    Dbxref A database cross reference.  See the section
	   "Ontology Associations and Db Cross References" for
	   details on the format.

    Ontology_term  A cross reference to an ontology term.  See
           the section "Ontology Associations and Db Cross References"
	   for details.

Multiple attributes of the same type are indicated by separating the
values with the comma "," character, as in:

       Parent=AF2312,AB2812,abc-3

Note that attribute names are case sensitive.  "Parent" is not the
same as "parent".

All attributes that begin with an uppercase letter are reserved for
later use.  Attributes that begin with a lowercase letter can be used
freely by applications.



=back

=cut

    ;



sub to_GFF3_format {
    my ($gene_obj, %preferences) = @_;
    
    my $gene_id = $gene_obj->{TU_feat_name};
    if ($gene_id =~ /;/) {
        $gene_id = "\"$gene_id\"";
    }
    
    my $strand = $gene_obj->get_orientation();
    
    my @noteText;
    
    if ($gene_obj->{is_pseudogene}) {
        push (@noteText, "(pseudogene)");
    }
    
    ## parse preferences
    my $asmbl_id = $preferences{seqid} || $gene_obj->{asmbl_id};
    my $source = $preferences{source} || $gene_obj->{source} || ".";
    
    unless (defined $asmbl_id) {
        confess "Error, no asmbl_id from gene_obj\n";
    }
    
    
    my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj->get_gene_span();
    my $com_name = $gene_obj->{com_name};
	unless ($com_name =~ /\w/) {
		$com_name = "";
	}
	
	if ($com_name) {
        if ($preferences{uri_encode_name}) {
            # uri escape it:
            use URI::Escape;
            $com_name = uri_escape($com_name);
        }
        else {
            unless (substr($com_name,0,1) =~ /\'|\"/ && substr($com_name, -1, 1) =~ /\'|\"/) {
                $com_name = "\"$com_name\"";
            }
        }
    }
    
	my $gene_alias = "";
	if (my $pub_locus = $gene_obj->{pub_locus}) {
		$gene_alias = "Alias=$pub_locus;";
	}
       
    my $feat_type = ($gene_obj->{gene_type} eq "protein-coding") ? "gene" : $gene_obj->{gene_type};
    

    my $gff3_text = "$asmbl_id\t$source\t$feat_type\t$gene_lend\t$gene_rend\t.\t$strand\t.\tID=$gene_id;Name=$com_name;$gene_alias\n";  ## note, non-coding gene features are currently represented by a simple single coordinate pair.
    
    if ($gene_obj->{gene_type} eq "protein-coding")  {
		
        my $gene_obj_ref = $gene_obj;
        
        foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms() ) {
            
            my $model_id = $gene_obj->{Model_feat_name};
            if ($model_id =~ /;/) {
                $model_id = "\"$model_id\"";
            }
            
            my $model_alias = "";
            if (my $model_locus = $gene_obj->{Model_pub_locus}) {
				$model_alias = "Alias=$model_locus;";
			}
			      
			my ($mrna_lend, $mrna_rend) = $gene_obj->get_transcript_span();
      
            $gff3_text .= "$asmbl_id\t$source\tmRNA\t$mrna_lend\t$mrna_rend\t.\t$strand\t.\tID=$model_id;Parent=$gene_id;Name=$com_name;$model_alias\n";
            
            ## mark the first and last CDS entries (for now, an unpleasant hack!)
            my @exons = $gene_obj->get_exons();
            ## find the first cds
            foreach my $exon (@exons) {
                if (my $cds = $exon->get_CDS_obj()) {
                    $cds->{first_cds} = 1;
                    last;
                }
            }
            @exons = reverse @exons;
            foreach my $exon (@exons) {
                if (my $cds = $exon->get_CDS_obj()) {
                    $cds->{last_cds} = 1;
                    last;
                }
            }
            
            my $prime5_partial = $gene_obj->is_5prime_partial();
            my $prime3_partial = $gene_obj->is_3prime_partial();
            
            
            ## annotate 5' utr
            if ($gene_obj->has_CDS() && $gene_obj->has_5prime_UTR()) {
                my @prime5_utr = $gene_obj->get_5prime_UTR_coords();
                if (@prime5_utr) {
                    my $utr_count = 0;
                    foreach my $coordset (@prime5_utr) {
                        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
                        $utr_count++;
                        my $utr_id = "$model_id.utr5p$utr_count";
                        $gff3_text .= "$asmbl_id\t$source\tfive_prime_UTR\t$lend\t$rend\t.\t$strand\t.\tID=$utr_id;Parent=$model_id\n";
                    }
                }
            }
            
			
			my $exon_counter = 0;
            foreach my $exon ($gene_obj->get_exons()) {
                $exon_counter++;
				my ($exon_lend, $exon_rend) = sort {$a<=>$b} $exon->get_coords();
                my $exon_ID_string = "";
                if (my $exon_feat_name = $exon->{feat_name}) {
                    $exon_ID_string = "$exon_feat_name";
                }
				else {
					$exon_ID_string = "$model_id.exon$exon_counter";
				}
                $gff3_text .= "$asmbl_id\t$source\texon\t$exon_lend\t$exon_rend\t.\t$strand\t.\tID=${exon_ID_string};Parent=$model_id\n";

                if (my $cds_obj = $exon->get_CDS_obj()) {
                    my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds_obj->get_coords();
                    my $phase = $cds_obj->{phase};
					if (defined($phase)) {
						## use GFF3 definition of phase, which is how many bases to trim before encountering first base of start
						if ($phase == 2) { 
							$phase = 1;
						}
						elsif ($phase == 1) {
							$phase = 2;
						}
						# phase 0 remains 0
					} 
					else {
						$phase =  "."; #use phase info if avail
                    }

					
					my $cds_ID_string = "cds.$model_id";
					
					# according to the GFF3 spec, CDS segments from the same coding region should have the same identifier.
					#if (my $cds_feat_name = $cds_obj->{feat_name}) {
					#	$cds_ID_string = "$cds_feat_name";
					#}
					#else {
					#	$cds_ID_string = "$model_id.cds$exon_counter";
					#}
					
                    my $partial_text = "";
                    if ($prime5_partial && $cds_obj->{first_cds}) {
                        $partial_text .= ";5_prime_partial=true";
                    }
                    if ($prime3_partial && $cds_obj->{last_cds}) {
                        $partial_text .= ";3_prime_partial=true";
                    }
                    
                    $gff3_text .= "$asmbl_id\t$source\tCDS\t$cds_lend\t$cds_rend\t.\t$strand\t$phase\tID=${cds_ID_string};Parent=$model_id$partial_text\n";
                }
            }
            
            ## annotate 3' utr
            if ($gene_obj->has_CDS() && $gene_obj->has_3prime_UTR()) {
                my @prime3_utr = $gene_obj->get_3prime_UTR_coords();
                if (@prime3_utr) {
                    my $utr_count = 0;
                    foreach my $coordset (@prime3_utr) {
                        my ($lend, $rend) = sort {$a<=>$b} @$coordset;
                        $utr_count++;
                        my $utr_id = "$model_id.utr3p$utr_count";
                        $gff3_text .= "$asmbl_id\t$source\tthree_prime_UTR\t$lend\t$rend\t.\t$strand\t.\tID=$utr_id;Parent=$model_id\n";
                    }
                }
                
            }
        }
        
    }  ## end of protein-coding genes
        

	## strip off any trailing whitespace and semicolons:
	my @lines = split (/\n/, $gff3_text);
	foreach my $line (@lines) {
		$line =~ s/\s+$//;
		$line =~ s/;$//;
	}
	
	$gff3_text = join ("\n", @lines) . "\n";
	
    return ($gff3_text);

}



=over 4

=item to_BED_format()

B<Description:> describes gene in BED format
B<Parameters:> (uri_encode => 1|0)
B<Returns:> string


BED format described here:
http://genome.ucsc.edu/FAQ/FAQformat.html#format1

	BED format


	

BED format provides a flexible way to define the data lines that are displayed in an annotation track. BED lines have three required fields and nine additional optional fields. The number of fields per line must be consistent throughout any single set of data in an annotation track. The order of the optional fields is binding: lower-numbered fields must always be populated if higher-numbered fields are used. 

The first three required BED fields are:

1. chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671). 

2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0. 

3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 

The 9 additional optional BED fields are:

4. name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode. 
	
5. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browsers translation of BED score values into shades of gray

6. strand - Defines the strand - either '+' or '-'. 

7. thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). 

8. thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays). 

9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser. 

10. blockCount - The number of blocks (exons) in the BED line. 

11. blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount. 

12. blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 

Example:
Heres an example of an annotation track that uses a complete BED definition:
track name=pairedReads description="Clone Paired Reads" useScore=1
chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601




=cut


sub to_BED_format {
	my $self = shift;
    my %params = @_;
    
	my $strand = $self->get_strand();

	my ($coding_lend, $coding_rend) = sort {$a<=>$b} $self->get_CDS_span();
	
	my $scaffold = $self->{asmbl_id};

	my $gene_id = $self->{TU_feat_name};
	my $trans_id = $self->{Model_feat_name};
	
	my $com_name = $self->{com_name} || "";

    my $score = $params{score} || 0;
    
	if (my $alias = $self->{pub_locus}) {
		$com_name = "Alias=$alias;$com_name";
	}
	
	
    if ($gene_id) {
        $com_name = "$gene_id;$com_name";
    }
    
	if ($trans_id) {
		$com_name = "ID=$trans_id;$com_name";
	}
	else {
        $com_name = "ID=$com_name";
    }
    
    if ($params{uri_encode}) {
        $com_name = uri_escape($com_name);
    }
    

	my @exons = sort {$a->{end5}<=>$b->{end5}} $self->get_exons();

    my @exon_coords;
	foreach my $exon (@exons) {
        
		my ($exon_lend, $exon_rend) = sort {$a<=>$b} $exon->get_coords();
        push (@exon_coords, [$exon_lend, $exon_rend]);
    }
    
	
	my @starts;
	my @lengths;

    my $gene_lend = $exon_coords[0]->[0];
    my $gene_rend = $exon_coords[$#exon_coords]->[1];
    
    foreach my $exon_coordset (@exon_coords) {
        my ($exon_lend, $exon_rend) = @$exon_coordset;
        
		my $start = $exon_lend - $gene_lend;
		push (@starts, $start);

		my $length = $exon_rend - $exon_lend + 1;
		push (@lengths, $length);
	}
		
	
	## construct bed output.

	$com_name =~ s/ /_/g;
	
	my $bed_line = join("\t", $scaffold, 
						$gene_lend-1, $gene_rend, 
						$com_name, 
						$score, 
						$strand,
						$coding_lend-1, $coding_rend,
						"0", # rgb info - use '.' to allow user customization in IGV.   Need 0 for compatibility with UCSC browser.
						scalar(@lengths),
						join(",", @lengths),
						join(",", @starts)
						) . "\n";

    foreach my $isoform ($self->get_additional_isoforms()) {
        $bed_line .= $isoform->to_BED_format(%params);
    }
    
	return($bed_line);
}



# static method, returns gene object.
sub BED_line_to_gene_obj {
	my ($bed_line) = @_;

	if (ref $bed_line) {
		confess "Error, static method, just provide bed text line, returns gene_obj";
	}
	
	
	my @x = split(/\t/, $bed_line);
	
	my $scaff = $x[0];
	my $gene_lend = $x[1] + 1;
	my $gene_rend = $x[2];

	my $com_name = $x[3];

	my $score = $x[4];
	my $orient = $x[5];
	
	if ($orient eq '*') {
		$orient = '+';
	}
	

	my $coding_lend = $x[6] + 1;
	my $coding_rend = $x[7];

	my $rgb_color = $x[8];

	my $num_exons = $x[9];

	my $lengths_text = $x[10];
	my $exon_relative_starts_text = $x[11];

	my @lengths = split(/,/, $lengths_text);
	my @exon_relative_starts = split(/,/, $exon_relative_starts_text);

	my @exons;

	while (@lengths) {
		my $len = shift @lengths;
		my $start = shift @exon_relative_starts;
		
		my $exon_lend = $gene_lend + $start;
		my $exon_rend = $exon_lend + $len - 1;
		

		print "Len: $len, start=$start   ====>  $exon_lend - $exon_rend\n" if $DEBUG;

		push (@exons, [$exon_lend, $exon_rend]);
		
	}

	
	print "Coding: $coding_lend-$coding_rend, Exons: " . Dumper (\@exons) if $DEBUG;
	
	my $gene_obj = new Gene_obj();
	$gene_obj->build_gene_obj_exons_n_cds_range(\@exons, $coding_lend, $coding_rend, $orient);
	
	$gene_obj->{com_name} = $com_name;
	$gene_obj->{asmbl_id} = $scaff;
	
	$com_name =~ s/\s+/\|/g; # reformat as an identifier with no whitespace
	
	$gene_obj->{TU_feat_name} = "$com_name";
	$gene_obj->{Model_feat_name} = "m.$com_name";
	
	return($gene_obj);
	
		
}






## Private, remove leading and trailing whitespace characters:
sub trim_leading_trailing_ws {
    my ($ref) = @_;
    if (ref $ref eq "SCALAR") {
        $$ref =~ s/^\s+|\s+$//g;
    } elsif (ref $ref eq "ARRAY") {
        foreach my $element (@$ref) {
            $element =~ s/^\s+|\s+$//g;
        }
    } else {
        my $type = ref $ref;
        die "Currently don't support trim_leading_trailing_ws(ref type: $type)\n";
    }
}



=over 4

=item to_GTF2_format()

B<Description:> provides gene in GTF2 format

B<Parameters:> genome_seq_ref, [properties_href]

B<Returns:> text


properties_href encodes preferences like so

    properties_href = { 
                          seqname => tigr_asmbl_id_1000, # by default, asmbl_id is used as encoded in gene_obj
                          
                          source => MyGenePrediction,   # by default, set to "TIGR"

                          include_comments => 0,   # turned on by default, indicating partial or pseudogenes with preceding comment lines

                      }



The GTF2 format is described here:
http://genes.cs.wustl.edu/GTF2.html

as follows:

GTF2 format (Revised Ensembl GTF)
Gene transfer format. This borrows from GFF, but has additional structure that warrants a separate definition and format name.
NEW! Validating Parser for GTF

Structure is as GFF, so the fields are:
<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

Here is a simple example with 3 translated exons. Order of rows is not important.

AB000381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";
AB000381 Twinscan  CDS          501   650   .   +   2  gene_id "001"; transcript_id "001.1";
AB000381 Twinscan  CDS          700   707   .   +   2  gene_id "001"; transcript_id "001.1";
AB000381 Twinscan  start_codon  380   382   .   +   0  gene_id "001"; transcript_id "001.1";
AB000381 Twinscan  stop_codon   708   710   .   +   0  gene_id "001"; transcript_id "001.1";

The whitespace in this example is provided only for readability. In GTF, fields must be separated by a single TAB and no white space.

<seqname>
The FPC contig ID from the Golden Path.

<source>
The source column should be a unique label indicating where the annotations came from --- typically the name of either a prediction program or a public database.

<feature>
The following feature types are required: "CDS", "start_codon", "stop_codon". The feature "exon" is optional, since this project will not evaluate predicted splice sites outside of protein coding regions. All other features will be ignored.

CDS represents the coding sequence starting with the first translated codon and proceeding to the last translated codon. Unlike Genbank annotation, the stop codon is not included in the CDS for the terminal exon.

<start> <end>
Integer start and end coordinates of the feature relative to the beginning of the sequence named in <seqname>.  <start> must be less than or equal to <end>. Sequence numbering starts at 1. Values of <start> and <end> that extend outside the reference sequence are technically acceptable, but they are discouraged for purposes of this project.

<score>
The score field will not be used for this project, so you can either provide a meaningful float or replace it by a dot.

<frame>
0 indicates that the first whole codon of the reading frame is located at 5'-most base. 1 means that there is one extra base before the first codon and 2 means that there are two extra bases before the first codon. Note that the frame is not the length of the CDS mod 3.

Here are the details excised from the GFF spec. Important: Note comment on reverse strand.

    '0' indicates that the specified region is in frame, i.e. that its first base corresponds to the first base of a codon. '1' indicates that there is one extra base, i.e. that the second base of the region corresponds to the first base of a codon, and '2' means that the third base of the region is the first base of a codon. If the strand is '-', then the first base of the region is value of <end>, because the corresponding coding region will run from <end> to <start> on the reverse strand.

[attributes]
All four features have the same two mandatory attributes at the end of the record:

    * gene_id value;     A globally unique identifier for the genomic source of the transcript
    * transcript_id value;     A globally unique identifier for the predicted transcript.

These attributes are designed for handling multiple transcripts from the same genomic region. Any other attributes or comments must appear after these two and will be ignored.

Attributes must end in a semicolon which must then be separated from the start of any subsequent attribute by exactly one space character (NOT a tab character).

Textual attributes should be surrounded by doublequotes.

Here is an example of a gene on the negative strand. Larger coordinates are 5' of smaller coordinates. Thus, the start codon is 3 bp with largest coordinates among all those bp that fall within the CDS regions. Similarly, the stop codon is the 3 bp with coordinates just less than the smallest coordinates within the CDS regions.

AB000123    Twinscan     CDS    193817    194022    .    -    2    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     CDS    199645    199752    .    -    2    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     CDS    200369    200508    .    -    1    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     CDS    215991    216028    .    -    0    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     start_codon   216026    216028    .    -    .    gene_id    "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     stop_codon    193814    193816    .    -    .    gene_id    "AB000123.1"; transcript_id "AB00123.1.2";

Note the frames of the coding exons. For example:

   1. The first CDS (from 216028 to 215991) always has frame zero.
   2. Frame of the 1st CDS =0, length =38.  (frame - length) % 3  = 1, the frame of the 2nd CDS.
   3. Frame of the 2nd CDS=1, length=140. (frame - length) % 3  = 2, the frame of the 3rd CDS.
   4. Frame of the 3rd CDS=2, length=108. (frame - length) % 3  =  2, the frame of the terminal CDS.
   5. Alternatively, the frame of terminal CDS can be calculated without the rest of the gene. Length of the terminal CDS=206. length % 3 =2, the frame of the terminal CDS.

Here is an example in which the "exon" feature is used. It is a 5 exon gene with 3 translated exons.

AB000381 Twinscan  exon         150   200   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  exon         300   401   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  CDS          380   401   .   +   0  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  exon         501   650   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  CDS          501   650   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  exon         700   800   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  CDS          700   707   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  exon         900  1000   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  start_codon  380   382   .   +   0  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
AB000381 Twinscan  stop_codon   708   710   .   +   0  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
  




=back

=cut



sub to_GTF2_format () {
    my $self = shift;
    my $genomic_seq_ref = shift;
    
    my $properties_href = shift;
    unless ($properties_href) {
        $properties_href = {};
    }
    

    ## need to adjust my frame definition so it's consistent with requirements above in spec.
    my $frame_convert = sub { 
        my $phase = shift;
        
        my %frame = ( 0 => 0,
                      1 => 2,
                      2 => 1 );
        return ($frame{$phase});
    };
    
    
    my $gtf2_text = "";
    
    my $gene_obj = $self;

    my $asmbl_id = $properties_href->{seqname} || $gene_obj->{asmbl_id} || die "Error, no asmbl_id as gene_obj att";
    
    my $source = $properties_href->{source} || "TIGR";
    
    my $gene_id = $gene_obj->{TU_feat_name};
    my $model_id = $gene_obj->{Model_feat_name};
    my $strand = $gene_obj->get_orientation();
    

    my $comment_line = "";
    if ($gene_obj->is_pseudogene()) {
        $comment_line .= "$model_id=pseudogene ";
    }
    
    if ( $gene_obj->{gene_type} eq "protein-coding") {
        
        if (! $gene_obj->is_pseudogene()) { 
            
            $gene_obj->set_CDS_phases($genomic_seq_ref);
            # also resets the 5' and 3' partiality attributes based on the longest orf.
            
            
            if ($gene_obj->is_5prime_partial()) {
                $comment_line .= "$model_id=5'partial ";
            }
            else {
                $gene_obj->validate_start_codon();
            }
            
            if ($gene_obj->is_3prime_partial() ) {
                $comment_line .= "$model_id=3'partial ";
            }      
            else {
                $gene_obj->validate_stop_codon();
            }
        }
        
        my @stop_codon_objs;
        my @start_codons;

        if (! $gene_obj->is_pseudogene()) {
            
            if (! $gene_obj->is_3prime_partial())  {
                @stop_codon_objs = $gene_obj->_remove_stop_codons();
                
                unless (@stop_codon_objs) {
                    confess $gene_obj->toString() . "Error, no stop codon objs retrieved for non 3' partial gene";
                }
            }
            if (! $gene_obj->is_5prime_partial()) {
                @start_codons = $self->_extract_start_codons();
                
                unless (@start_codons) {
                    confess $gene_obj->toString() . "Error, no start codon extracted for non 5'partial gene.";
                }
            }
            
        }
        
        foreach my $start_codon (@start_codons) {
            my ($start_lend, $start_rend) = sort {$a<=>$b} $start_codon->get_coords();
            my $phase = &$frame_convert($start_codon->{phase});
            $gtf2_text .= "$asmbl_id\t$source\tstart_codon\t$start_lend\t$start_rend\t.\t$strand\t$phase\tgene_id \"$gene_id\"; transcript_id \"$model_id\";\n";
        }
        
        
        foreach my $exon ($gene_obj->get_exons()) {
            my ($exon_lend, $exon_rend) = sort {$a<=>$b} $exon->get_coords();
            $gtf2_text .= "$asmbl_id\t$source\texon\t$exon_lend\t$exon_rend\t.\t$strand\t.\tgene_id \"$gene_id\"; transcript_id \"$model_id\";\n";
            
            if ($gene_obj->is_pseudogene()) { next; } # don't bother trying to report nonsensical CDSs.
            
            if (my $cds_obj = $exon->get_CDS_obj()) {
                my ($cds_lend, $cds_rend) = sort {$a<=>$b} $cds_obj->get_coords();
                my $phase = $cds_obj->{phase};
                unless (defined($phase)) {
                    die "Error, no phase defined for cds($cds_lend-$cds_rend) of gene" . $gene_obj->toString();
                }
                $phase = &$frame_convert($phase);

                $gtf2_text .= "$asmbl_id\t$source\tCDS\t$cds_lend\t$cds_rend\t.\t$strand\t$phase\tgene_id \"$gene_id\"; transcript_id \"$model_id\";\n";
            }
        }
        
        foreach my $stop_codon (@stop_codon_objs) {
            my ($stop_lend, $stop_rend) = sort {$a<=>$b} $stop_codon->get_coords();
            my $phase = &$frame_convert($stop_codon->{phase});
            $gtf2_text .= "$asmbl_id\t$source\tstop_codon\t$stop_lend\t$stop_rend\t.\t$strand\t$phase\tgene_id \"$gene_id\"; transcript_id \"$model_id\";\n";
        }
        
        foreach my $isoform ($gene_obj->get_additional_isoforms() ) {
            $gtf2_text .= $isoform->to_GTF2_format($genomic_seq_ref, $properties_href);
        }
    }

    if ($comment_line) {
        # prefix with \# to actually comment it in the file
        $comment_line = "#$comment_line\n";
    }
    
    my $comment_flag = $properties_href->{include_comments};
    if (defined ($comment_flag) && $comment_flag == 0) {
        $comment_line = ""; # clear it
    }
    
    
    return ($comment_line . $gtf2_text);
}


sub _extract_start_codons {
    my $self = shift;

    ## 5' partiality attribute is trusted here !!!
    
    if ($self->is_5prime_partial()) {
        return();
    }

    my @exons = $self->get_exons();
    my $orientation = $self->get_orientation();
    
    my @start_codons;
    
    my $found_cds_flag = 0;
    
    for (my $i = 0; $i <= $#exons; $i++) {
        if (my $cds = $exons[$i]->get_CDS_obj()) {
            # found first cds
            $found_cds_flag = 1;
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            my $cds_len = $cds->length();
            if ($cds_len >= 3) {
                ## got start codon in entirety
                if ($orientation eq '+') {
                    push (@start_codons, CDS_exon_obj->new($cds_end5, $cds_end5+2)->set_phase(0));
                    last;
                } 
                else {
                    push (@start_codons, CDS_exon_obj->new($cds_end5, $cds_end5-2)->set_phase(0));
                    last;
                }
            }
            else {
                ## split start codon
                push (@start_codons, $cds); # add current cds as start codon part
                my $missing_length = 3 - $cds_len;
                
                ## examine next cds exon for part of it:
                my $next_cds = $exons[$i+1]->get_CDS_obj();
                unless (ref $next_cds) {
                    die "Error, no next cds for split start codon" . $self->toString();
                }
                
                my ($next_cds_end5, $next_cds_end3) = $next_cds->get_coords();
                my $next_cds_len = $next_cds->length();
                
                if ($next_cds_len >= $missing_length) {
                    # great, this has everything we need
                    if ($orientation eq '+') {
                        push (@start_codons, 
                              CDS_exon_obj->new($next_cds_end5, $next_cds_end5 + $missing_length-1)->set_phase($next_cds->{phase}));
                        last;
                    }
                    else {
                        push (@start_codons, 
                              CDS_exon_obj->new($next_cds_end5, $next_cds_end5 - $missing_length + 1)->set_phase($next_cds->{phase}));
                        last;
                    }
                }
                else {
                    ## another split start codon portion.  Just add the current cds, and get the first bp from the next cds
                    push (@start_codons, $next_cds);
                    
                    my $final_cds = $exons[$i+2]->get_CDS_obj();
                    unless (ref $final_cds) {
                        die "Error getting final cds of three-part split start codon";
                    }
                    unless ($final_cds->{phase} == 2) {
                        die "Error, final cds of three-part stop codon is not in phase 2 ";
                    }
                    my ($final_cds_end5, $final_cds_end3) = $final_cds->get_coords();
                    push (@start_codons,
                          CDS_exon_obj->new($final_cds_end5, $final_cds_end5)->set_phase(2));
                    last;
                }
            } # end of split start codon
            
        } # end of found cds
    } # end of foreach exon
    
    unless ($found_cds_flag) {
        die "Error, no cds exon found in search of start codon";
    }
    
    unless (@start_codons) {
        die "Error, no start codons found";
    }
    ## ensure start codons sum to 3
    my $sum_len = 0;
    foreach my $start_codon (@start_codons) {
        $sum_len += $start_codon->length();
    }
    unless ($sum_len == 3) {
        print "Error, sum len of start codons != 3 ( = $sum_len, instead) " . $self->toString() . "starts:\n";
        my $i=0;
        foreach my $start (@start_codons) {
            $i++;
            print "start($i): " . $start->toString();
        }
        die;
    }
    
    return (@start_codons);
    
}



sub _remove_stop_codons {
    my $self = shift;
    
    ## 3' partiality attribute is trusted here !!!
    
    if ($self->is_3prime_partial()) {
        return ();
    }
    
    my $orientation = $self->get_orientation();
    my @exons = reverse $self->get_exons(); # examining exons in reverse order, starting from stop codon direction.
    
    my @stop_codons;
    
    my $found_cds_flag = 0;

    ## find first exon
    for (my $i=0; $i <= $#exons; $i++) {
        if (my $cds = $exons[$i]->get_CDS_obj()) {
            
            $found_cds_flag = 1;
            
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            
            my $cds_length = $cds->length();
            if ($cds_length > 3) {
                ## cds exon encodes more than just the stop codon
                if ($orientation eq '+') {
                    $cds->{end3} -= 3;
                    push (@stop_codons, CDS_exon_obj->new($cds_end3 - 2, $cds_end3)->set_phase(0));
                }
                else {
                    $cds->{end3} += 3;
                    push (@stop_codons, CDS_exon_obj->new($cds_end3 + 2, $cds_end3)->set_phase(0));
                }
                last;
                
            }
            elsif ($cds_length == 3) {
                ## Just a stop codon exon.  We can remove it.
                push (@stop_codons, $cds);
                $exons[$i]->{CDS_exon_obj} = 0; # nullified
                last;
            }
            
            else {
                ## cds exon encodes a split stop codon
                push (@stop_codons, $cds); # just add the last portion of stop codon
                $exons[$i]->{CDS_exon_obj} = 0; # nullified
                
                ## check next portion of cds exon to see if it contains the rest of the stop 
                my $next_exon = $exons[$i+1];
                unless (ref $next_exon) {
                    die "Error, incomplete stop codon and not enough exons! ";
                }
                my $missing_stop_length = 3 - $cds_length;
                my $next_cds_obj = $next_exon->get_CDS_obj();
                unless (ref $next_cds_obj) {
                    die "Error, next cds obj is missing!";
                }
                
                my $next_cds_length = $next_cds_obj->length();
                my ($cds_end5, $cds_end3) = $next_cds_obj->get_coords();
                if ($next_cds_length <= $missing_stop_length) {
                    ## encodes only the second part of the stop codon
                    # add and nullify
                    push (@stop_codons, $next_cds_obj);
                    $next_exon->{CDS_exon_obj} = 0;
                    
                    ## get the very last part of the stop 
                    $missing_stop_length -= $next_cds_length;
                    if ($missing_stop_length > 0) {
                        ## must be still missing the first bp of the stop codon
                        if ($missing_stop_length != 1) {
                            die "Error, too much of the stop codon is left (missing_length = $missing_stop_length).  Should only be 1 ";
                        }
                        my $next_exon = $exons[$i+2];
                        unless (ref $next_exon) {
                            die "Error, second next exon is unavail ";
                        }
                        my $next_cds_obj = $next_exon->get_CDS_obj();
                        unless (ref $next_cds_obj) {
                            die "Error, second next cds obj is unavail";
                        }
                        my $cds_length = $next_cds_obj->length();
                        my ($cds_end5, $cds_end3) = $next_cds_obj->get_coords();
                        if ($cds_length > 1) {
                            if ($orientation eq '+') {
                                $next_cds_obj->{end3}-=1;
                                push (@stop_codons, CDS_exon_obj->new($cds_end3, $cds_end3)->set_phase(0));
                            }
                            else {
                                $next_cds_obj->{end3}+=1;
                                push (@stop_codons, CDS_exon_obj->new($cds_end3, $cds_end3)->set_phase(0));
                            }
                        }
                    } 
                }
                else {
                    # split stop codon
                    #missing length of cds exon is present in the second portion of the stop 
                    if ($orientation eq '+') {
                        $next_cds_obj->{end3} -= $missing_stop_length;
                        push (@stop_codons, CDS_exon_obj->new($cds_end3 - $missing_stop_length + 1, $cds_end3)->set_phase(0));
                    }
                    else {
                        $next_cds_obj->{end3} += $missing_stop_length;
                        push (@stop_codons, CDS_exon_obj->new($cds_end3 + $missing_stop_length -1, $cds_end3)->set_phase(0));
                    }
                }
            } # end of split stop codon
            
            last;
            
        } # end of found cds obj
        

    } # end of foreach exon 
    
    unless ($found_cds_flag) {
        die "Error, no cds exon was found. ";
    }


    unless (@stop_codons) {
        die "Error, no stop codons extracted from non-partial gene.";
    }

    @stop_codons = reverse @stop_codons; # reorder according to gene direction

    ## make sure sum (stop_codons) length == 3
    my $sum_len = 0;
    foreach my $stop_codon (@stop_codons) {
        $sum_len += $stop_codon->length();
    }
    if ($sum_len != 3) {
        print "Error, stop codons sum length != 3 ( = $sum_len, instead) " . $self->toString();
        my $i=0;
        foreach my $stop_codon (@stop_codons) {
            print "stop($i): " . $stop_codon->toString();
        }

        die;
    }
    
    return (@stop_codons);
    
}
                


sub has_CDS {
    my $self = shift;

    foreach my $exon ($self->get_exons()) {
        if (ref ($exon->get_CDS_obj())) {
            return (1);
        }
    }

    return (0); # no cds entry found
}


####
sub set_CDS_phases_from_init_phase {
    my ($self, $init_phase) = @_;

    my @exons = $self->get_exons();

    my $curr_cds_len = $init_phase;
    
    foreach my $exon (@exons) {
        if (my $cds = $exon->get_CDS_obj()) {
            $cds->set_phase($curr_cds_len % 3);
            my $cds_len = $cds->length();
            $curr_cds_len += $cds_len;
        }
    }

    return;
}



sub set_CDS_phases {
    my ($self, $genomic_seq_ref) = @_;
        

    my $start_pos = 1;
    if ($self->has_CDS() && ! $self->is_pseudogene()) {
     
        $self->create_all_sequence_types($genomic_seq_ref);
        
        my $cds_sequence = $self->get_CDS_sequence();
        my $protein_seq = $self->get_protein_sequence();
        
        ## first, clear the partial attributes:
        $self->set_5prime_partial(0);
        $self->set_3prime_partial(0);
        
        
        if ($protein_seq !~ /^M/) {
            # lacks start codon
            $self->set_5prime_partial(1);
        }
        if ($protein_seq !~ /\*$/) {
            # lacks stop codon
            $self->set_3prime_partial(1);
        }
        
        $start_pos = $self->_get_cds_start_pos($cds_sequence);
        
        ## must set phase based on codon start position:
        ## my definition of phase here is the actual codon position of the first base in the CDS sequence.
        ##  (note this differs from the GFF3 spec, and is adjusted for in the to_GFF3_format() method.
        ##
                
        my $first_phase;
        if ($start_pos == 0) {
            # ATG XXX ...
            # 012 012 012
            
            $first_phase = 0;  
        }
        elsif ($start_pos == 1) {
            # XAT GXX ...
            # 201 201 201
            
            $first_phase = 2;
        }
        elsif ($start_pos == 2) {
            # XXA TGX ...
            # 120 120 120 

            $first_phase = 1
        }
        else {
            confess "Error, start pos: $start_pos doesn't make sense here... must be a bug.";
        }

        my @exons = $self->get_exons();
        my @cds_objs;
        foreach my $exon (@exons) {
            my $cds = $exon->get_CDS_obj();
            if (ref $cds) {
                push (@cds_objs, $cds);
            }
        }
        
        my $cds_obj = shift @cds_objs;
        $cds_obj->{phase} = $first_phase;
        my $cds_length = abs ($cds_obj->{end3} - $cds_obj->{end5}) + 1;

        if ($first_phase != 0) {
            # and now I understand why the GFF3 phase definition differs from mine. :-)  
            if ($first_phase == 1) {
                $cds_length -= 2;
            }
            elsif ($first_phase == 2) {
                $cds_length -= 1;
            }
            else {
                confess "Error, first phase set to: $first_phase, which is nonsensical";
            }
        }
                
        while (@cds_objs) {
            my $next_cds_obj = shift @cds_objs;
            $next_cds_obj->{phase} = $cds_length % 3;
            $cds_length += abs ($next_cds_obj->{end3} - $next_cds_obj->{end5}) + 1;
        }
    }
    
    foreach my $isoform ($self->get_additional_isoforms()) {
        $isoform->set_CDS_phases($genomic_seq_ref);
    }
    
    return;

}


sub get_first_CDS_segment {
    my $gene_obj = shift;
    my @exons = $gene_obj->get_exons();
    
    foreach my $exon (@exons) {
        if (my $cds = $exon->get_CDS_exon_obj()) {
            return ($cds);
        }
    }

    return undef;
}

sub _get_cds_start_pos {
    my ($self, $cds_sequence) = @_;
    my $cds_length = length($cds_sequence);
    # if cds is set of triplets, assume translate at codon pos 1.
    my $codon_start;
    
    ## must determine where translation starts:
    my $new_orfFinder = new Longest_orf();
    $new_orfFinder->allow_partials();
    $new_orfFinder->forward_strand_only();
    
	my $longest_orf = $new_orfFinder->get_longest_orf($cds_sequence);
   			
	unless (ref $longest_orf) {
		die "No longest ORF found in sequence";
	}

	## examine the first three ORFs, prefer long orf with stop codon.
    my $orfPos = $longest_orf->{start}; #init to first, longest orf.
	unless (defined $orfPos) {
		die "Error, orfPos not defined! " . Dumper ($longest_orf);
	}
	
	my $bestOrfPos;
    my @allOrfs = $new_orfFinder->orfs();
	    
    for my $orfIndex (0..2) {
        my $orf = $allOrfs[$orfIndex];
        if ($orf) {
            my $start = $orf->{start};
            my $length = $orf->{length};
            my $protein = $orf->{protein};
            if ($length > $cds_length - 3 && $start <= 3 && $protein =~ /\*$/) {
                unless ($bestOrfPos) {
                    $bestOrfPos = $start;
                }
            }
        }
    }
    
    if ($bestOrfPos && $bestOrfPos != $orfPos) {
        $orfPos = $bestOrfPos;
    }
    
    if ($orfPos >3) {
        confess "Error, longest ORF is found at position $orfPos, and should be between 1 and 3.  What's wrong with your gene?" . $self->toString();
    }
    $codon_start = $orfPos;

    #longest orf apparently using 1-based rather than 0-based coordinates.
        
    $codon_start -= 1;  
    
    return ($codon_start);
}


=over 4
        
=item dispose()

B<Description:> Sets all attributes = 0, hopefully to faciliate targeting for garbage collection. (experimental method) 

B<Parameters:> none

B<Returns:> none

=back

=cut

sub dispose {
    my $self = shift;
    foreach my $att (keys %$self) {
	$self->{$att} = 0;
    }
}



sub DESTROY {
    my $self = shift;

    warn "DESTROYING gene_obj: " . $self->{TU_feat_name} . "," . $self->{Model_feat_name} . "\n" if $main::DEBUG;

}


sub validate_start_codon {
    ## requires that you have the CDS sequence already set
    my $self = shift;

    my $cds_sequence = $self->get_CDS_sequence() or confess "Error, cannot get CDS sequence.  It must be built prior to calling this method";
    ## currently, only trust Met start codons.
    my $start_codon = uc substr($cds_sequence, 0, 3);
    if ($start_codon ne "ATG") {
        die $self->toString() . "Error, start codon is not M (codon $start_codon instead)!";
        # call within an eval block to catch exception
    }
}


sub validate_stop_codon {
    ## requires that you have the CDS sequence already set
    my $self = shift;

    my $cds_sequence = $self->get_CDS_sequence() or confess "Error, cannot get CDS sequence.  It must be built prior to calling this method";
    
    my @stop_codons = &Nuc_translator::get_stop_codons();
    
    my $curr_stop_codon = substr($cds_sequence, length($cds_sequence)-3, 3);
    
    my $found_stop_codon_flag = 0;
    foreach my $stop (@stop_codons) {
        if ($stop eq $curr_stop_codon) {
            $found_stop_codon_flag = 1;
            last;
        }
    }

    unless ($found_stop_codon_flag) {
        die $self->toString() . "Error, stop codon $curr_stop_codon is not an acceptable stop codon: [@stop_codons]\n";
    }

}




######################################################################################################################################
######################################################################################################################################


=head1 NAME

package mRNA_exon_obj

=cut

=head1 DESCRIPTION

    The mRNA_exon_obj represents an individual spliced mRNA exon of a gene.  The coordinates of the exon can be manipulated, and the mRNA_exon_obj can contain a single CDS_exon_obj.  A mRNA_exon_obj lacking a CDS_exon_obj component is an untranslated (UTR) exon.

    A mature Gene_obj is expected to have at least one mRNA_exon_obj component.

=cut


package mRNA_exon_obj;

use strict;
use warnings;
use Storable qw (store retrieve freeze thaw dclone);

=over 4

=item new()

B<Description:> Instantiates an mRNA_exon_obj

B<Parameters:> <(end5, end3)>

The end5 and end3 coordinates can be optionally passed into the constructor to set these attributes.  Alternatively, the set_coords() method can be used to set these values.

B<Returns:> $mRNA_exon_obj

=back

=cut


    ;

sub new {
    shift;
    my $self = { end5 => 0,   # stores end5 of mRNA exon
                 end3 => 0,   # stores end3 of mRNA exon
                 CDS_exon_obj => 0,   # stores object reference to CDS_obj
                 feat_name => 0,    # stores TIGR temp id
                 strand => undef,   #   +|-
                 };
    
    # end5 and end3 can be included as parameters in constructor.
    if (@_) {
        my ($end5, $end3) = @_;
        if (defined($end5) && defined($end3)) {
            $self->{end5} = $end5;
            $self->{end3} = $end3;
        }
    }
    
    bless ($self);
    return ($self);
}



=over 4

=item get_CDS_obj()

B<Description:> Retrieves the CDS_exon_obj component of this mRNA_exon_obj

B<Parameters:> none

B<Returns:> $cds_exon_obj

If no CDS_exon_obj is attached, returns 0

=back

=cut

    ;

sub get_CDS_obj {
    my $self = shift;
    return ($self->{CDS_exon_obj});
}


## alias
sub get_CDS_exon_obj {
    my $self = shift;
    return ($self->get_CDS_obj());
}


=over 4

=item get_mRNA_exon_end5_end3()

B<Description:> Retrieves the end5, end3 coordinates of the exon

**Method Deprecated**, use get_coords()

B<Parameters:> none

B<Returns:> (end5, end3)

=back

=cut


sub get_mRNA_exon_end5_end3 {
    my $self = shift;
    return ($self->{end5}, $self->{end3});
}



=over 4

=item set_CDS_exon_obj()

B<Description:> Sets the CDS_exon_obj of the mRNA_exon_obj

B<Parameters:> $cds_exon_obj

B<Returns:> none

=back

=cut

    ;
sub set_CDS_exon_obj {
    my $self = shift;
    my $ref = shift;
    if (ref($ref)) {
        $self->{CDS_exon_obj} = $ref;
    }
}



####
sub delete_CDS_exon_obj {
    my $self = shift;
    $self->{CDS_exon_obj} = undef;
    return;
}


=over 4

=item add_CDS_exon_obj()

B<Description:> Instantiates and adds a new CDS_exon_obj to the mRNA_exon_obj given the CDS coordinates.

B<Parameters:> (end5, end3)

B<Returns:> none

=back

=cut


sub add_CDS_exon_obj {
    my $self = shift;
    my ($end5, $end3) = @_;
    my $cds_obj = CDS_exon_obj->new ($end5, $end3);
    $self->set_CDS_exon_obj($cds_obj);
}


=over 4

=item set_feat_name()

B<Description:> Sets the feat_name attribute of the mRNA_exon_obj

B<Parameters:> $feat_name

B<Returns:> none

=back

=cut



sub set_feat_name {
    my $self = shift;
    my $feat_name = shift;
    $self->{feat_name} = $feat_name;
}


=over 4

=item clone_exon()

B<Description:> Creates a deep clone of this mRNA_exon_obj, using dclone() of Storable.pm

B<Parameters:> none

B<Returns:> $mRNA_exon_obj

=back

=cut
    
    

sub clone_exon {
    my $self = shift;
  
    my $clone_exon = dclone($self);
        
    return ($clone_exon);
}



=over 4

=item get_CDS_end5_end3 ()

B<Description:> Retrieves end5, end3 of the CDS_exon_obj component of this mRNA_exon_obj

B<Parameters:> none

B<Returns:> (end5, end3)

An empty array is returned if no CDS_exon_obj is attached.

=back

=cut


sub get_CDS_end5_end3 {
    my $self = shift;
    my $cds_obj = $self->get_CDS_obj();
    if ($cds_obj) {
        return ($cds_obj->get_CDS_end5_end3());
    } else {
        return ( () );
    }
}


=over 4

=item get_coords()

B<Description:> Retrieves the end5, end3 coordinates of this mRNA_exon_obj

B<Parameters:> none

B<Returns:> (end5, end3)

=back

=cut


sub get_coords {
    my $self = shift;
    return ($self->get_mRNA_exon_end5_end3());
}


=over 4

=item set_coords()

B<Description:> Sets the end5, end3 coordinates of the mRNA_exon_obj

B<Parameters:> (end5, end3)

B<Returns:> none

=back

=cut


## simpler coord setting (end5, end3)
sub set_coords {
    my $self = shift;
    my $end5 = shift;
    my $end3 = shift;
    $self->{end5} = $end5;
    $self->{end3} = $end3;
}


=over 4

=item get_strand()

B<Description:> Retrieves the orientation of the mRNA_exon_obj based on gene models transcribed orientation.

B<Parameters:> none

B<Returns:> +|-|undef

If end5 == end3, strand orientation cannot be inferred based on coordinates alone, so undef is returned.

=back

=cut


    ;

sub get_orientation {
    # determine positive or reverse orientation
    my $self = shift;
    return ($self->{strand});
}


sub get_strand { ## preferred
	my $self = shift;
	return($self->get_orientation());
}


####
sub merge_exon {
    my $self = shift;
    my $other_exon = shift;

    my $cds = $self->get_CDS_exon_obj();
    
    my $other_cds = $other_exon->get_CDS_exon_obj();

    if ($other_cds) {
        if ($cds) {
            $cds->merge_CDS($other_cds);
        }
        else {
            # current exon lacks cds. Set this one to it.
            $self->set_CDS_exon_obj($other_cds);
        }
    }


    ## merge the exons.
    my @coords = sort {$a<=>$b} ($self->get_coords(), $other_exon->get_coords());
    my $lend = shift @coords;
    my $rend = pop @coords;
    
    my ($new_end5, $new_end3) = ($self->get_orientation() eq '+') ? ($lend, $rend) : ($rend, $lend);

    $self->set_coords($new_end5, $new_end3);

    return;
}
    





=over 4

=item toString()

B<Description:> Provides a textual description of the mRNA_exon_obj 

B<Parameters:> none

B<Returns:> $text

=back

=cut

    ;


sub toString {
    my $self = shift;
    my @coords = $self->get_mRNA_exon_end5_end3();
    my $feat_name = $self->{feat_name};
    my $text = "";
    if ($feat_name) {
        $text .= "feat_name: $feat_name\t";
    }
    $text .= "end5 " . $coords[0] . "\tend3 " . $coords[1] . "\n";
    return ($text);
}


sub length {
    my $self = shift;

    my $len = abs ($self->{end5} - $self->{end3}) + 1;
    
    return($len);
}





##########################################################################################################################
##########################################################################################################################



=head1 NAME

package CDS_exon_obj

=cut


=head1 DESCRIPTION

    The CDS_exon_obj represents the protein-coding portion of an mRNA_exon_obj.

=cut



package CDS_exon_obj;

use strict;
use warnings;
use Storable qw (store retrieve freeze thaw dclone);
use Carp;


=over 4

=item new()

B<Description:>  Cosntructor for the CDS_exon_obj

B<Parameters:> <(end5, end3)>

The (end5, end3) parameter is optional.  Alternatively, the set_coords() method can be used to set these values.

B<Returns:> $cds_exon_obj

=back

=cut

    ;

sub new {
    shift;
    my $self = { end5 => 0,   #stores end5 of cds exon
                 end3 => 0,    #stores end3 of cds exon
                 phase => undef, #must set if to output in gff3 format.
                 feat_name => 0, #tigr's temp id
                 strand => undef,   # +|-
             };
    
    
    # end5 and end3 are allowed constructor parameters
    if (@_) {
        my ($end5, $end3) = @_;
        if (defined ($end5) && defined ($end3)) {
            $self->{end5} = $end5;
            $self->{end3} = $end3;
        }
    }
    bless ($self);
    return ($self);
}



=over 4

=item set_feat_name()

B<Description:> Sets the feat_name attribute value of the CDS_exon_obj 

B<Parameters:> $feat_name

B<Returns:> none

=back

=cut


sub set_feat_name {
    my $self = shift;
    my $feat_name = shift;
    $self->{feat_name} = $feat_name;
}


=over 4

=item get_CDS_end5_end3()

B<Description:> Retrieves the end5, end3 coordinates of the CDS_exon_obj

** Method deprecated **, use get_coords()


B<Parameters:> none

B<Returns:> (end5, end3)

=back

=cut


sub get_CDS_end5_end3 {
    my $self = shift;
    return ($self->{end5}, $self->{end3});
}



=over 4

=item set_coords()
    
B<Description:> Sets the (end5, end3) values of the CDS_exon_obj 

B<Parameters:> (end5, end3)

B<Returns:> none

=back

=cut



sub set_coords {
    my $self = shift;
    my $end5 = shift;
    my $end3 = shift;
    $self->{end5} = $end5;
    $self->{end3} = $end3;
}

=over 4

=item get_coords()

B<Description:> Retrieves the (end5, end3) coordinates of the CDS_exon_obj

B<Parameters:> none

B<Returns:> (end5, end3)


The get_coords() method behaves similarly among Gene_obj, mRNA_exon_obj, and CDS_exon_obj, and is generally preferred to other existing methods for extracting these coordinate values.  Other methods persist for backwards compatibility with older applications, but have been largely deprecated.


=back

=cut



sub get_coords {
    my $self = shift;
    return ($self->get_CDS_end5_end3());
}


=over 4

=item get_orientation()

B<Description:> Retrieves the orientation of the CDS_exon_obj based on gene models orientation.

B<Parameters:> none

B<Returns:> +|-|undef

undef returned if end5 == end3

=back

=cut

    ;

sub get_orientation {
    # determine positive or reverse orientation
    my $self = shift;
    return ($self->{strand});
}


sub get_strand { ## preferred
	my $self = shift;
	return($self->get_orientation());
}


=over 4

=item toString()

B<Description:> Retrieves a textual description of the CDS_exon_obj

B<Parameters:> none

B<Returns:> $text

=back

=cut



=over 4

=item clone_cds()

B<Description:> Creates a deep clone of this CDS_exon_obj, using dclone() of Storable.pm

B<Parameters:> none

B<Returns:> $mRNA_exon_obj

=back

=cut
    
    

sub clone_cds {
    my $self = shift;
  
    my $clone_cds = dclone($self);
        
    return ($clone_cds);
}


=over 4

=item length()

B<Description:> length of this cds segment

B<Parameters:> none

B<Returns:> int

=back

=cut
    

sub length {
    my $self = shift;
    my $length = abs ($self->{end3} - $self->{end5}) + 1;
    return ($length);
}



=over 4

=item  set_phase()

B<Description:> set phase of the CDS incident bp

B<Parameters:> [012]

B<Returns:> self


phase 0 = first bp of codon
phase 1 = second bp of codon
phase 2 = third bp of codon


=back

=cut
    

sub set_phase {
    my $self = shift;
    my $phase = shift;
    $self->{phase} = $phase;
    return($self);
}

=over 4

=item  get_phase()

B<Description:> gets phase of the CDS incident bp

B<Parameters:> none

B<Returns:> [012] or undef if not set


phase 0 = first bp of codon
phase 1 = second bp of codon
phase 2 = third bp of codon


=back

=cut
    



sub get_phase {
    my $self = shift;
    my $phase = $self->{phase};
    return($phase);
}



####
sub merge_CDS {
    my $self = shift;
    my $other_cds = shift;
    
    my $orientation = $self->get_orientation();
    unless ($orientation) {
        confess "Error, self CDS lacks orientation\n";
    }
    
    my @coords = sort {$a<=>$b} ($self->get_coords(), $other_cds->get_coords());
    my $lend = shift @coords;
    my $rend = pop @coords;

    unless ($lend && $rend) {
        confess "Error, trying to merge CDSs but coordinates are not available: \n"
            . "self: " . $self->toString()
            . "\n"
            . "other: " . $other_cds->toString() . "\n";
    }
    
    my ($end5, $end3) = ($orientation eq '+') ? ($lend, $rend) : ($rend, $lend);
    
    $self->set_coords($end5, $end3);
}

sub toString {
    my $self = shift;
    my @coords = $self->get_CDS_end5_end3();
    my $feat_name = $self->{feat_name};
    my $text = "";
    if ($feat_name) {
        $text .= "feat_name: $feat_name\t";
    }
    $text .= "end5 " . $coords[0] . "\tend3 " . $coords[1] . "\n";
    return ($text);
}


1;














