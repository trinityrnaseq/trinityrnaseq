
###########################
### Class Alignment_segment
###########################

=head1 NAME

package CDNA::Alignment_segment

=head1 DESCRIPTION

Provides an object representation of alignment segments which are built into a single CDNA_alignment object.

=cut



package CDNA::Alignment_segment;
use strict;
use Data::Dumper;


=over 4

=item new()

B<Description:> Instantiates a new Alignment_segment object.

B<Parameters:> $genomic_end5, $genomic_end3, $cdna_end5, $cdna_end3, $per_id

B<Returns:> Alignment_segment_obj

Alignment_segment_obj is an object of type CDNA::Alignment_object

Use the methods described below.  In addition, the following fields are supported:

B<orientation>  (orientation of the alignment segment)

B<lend>  (left end of the alignment segment corresponding to the genomic sequence)

B<rend>  (right end of the alignmetn segment corresponding to the genomic sequence)

note: lend <= rend in all cases; must use the B<orientation> field to determine cDNA alignment orientation for the segment.

B<mlend> (the cDNA coordinate corresponding to the lend alignment coordinate)

B<mrend> (the cDNA coordinate corresponding to the rend alignment coordinate)

=back

=cut


sub new {
    my $packagename = shift;
    my ($genomic_end5, $genomic_end3, $cdna_end5, $cdna_end3, $per_id) = @_;
    my $orientation = '?'; #initialize
    my ($lend, $rend, $mlend, $mrend) = ($genomic_end5, $genomic_end3, $cdna_end5, $cdna_end3);
    
    #reorient coordsets so that cdna coordinates are always in forward orientation.
    if ($cdna_end5 > $cdna_end3) { #swap coordsets
        ($lend, $rend) = ($rend, $lend);
        ($mlend, $mrend) = ($mrend, $mlend);
    }
    
    ## Check orientation and adjust lend, rend accordingly.
    if ($lend > $rend) {
        $orientation = '-';
        ($lend, $rend) = ($rend, $lend);
        ($mlend, $mrend) = ($mrend, $mlend);
    } elsif ($lend < $rend) { #keep coords way they are.
        $orientation = '+';
    }
    
    my $self = {
        orientation=>$orientation, ## should be [+-]
        lend=>$lend,
        rend=>$rend,
        mlend=>$mlend,  ## cDNA coordinate that maps to lend of alignment.
        mrend=>$mrend,
        per_id => $per_id,
        type=>undef(), # [first|last|internal|single]
        has_left_splice_junction=>0, #flag indicating whether the consensus is present.
        has_right_splice_junction=>0,
        left_splice_site_chars=>undef(), #store the two characters at that splice junction.
        right_splice_site_chars=>undef()
        };
    bless ($self, $packagename);
    return ($self);
}

sub set_coords {
    my $self = shift;
    my ($c1, $c2) = @_;
    ($c1, $c2) = sort {$a<=>$b} ($c1, $c2);
    $self->{lend} = $c1;
    $self->{rend} = $c2;
}


####
sub get_aligned_orientation {
    my $self = shift;
    return ($self->{orientation});
}



=over 4

=item get_coords()

B<Description:> Retrieves the lend, rend for the alignment segment.

B<Parameters:> none.

B<Returns:> ($lend, $rend)

=back

=cut


sub get_coords {
    my $self = shift;
    return ($self->{lend}, $self->{rend});
}

# private.
sub set_mcoords () {
    my $self = shift;
    my ($mlend, $mrend) = @_;
    $self->{mlend} = $mlend;
    $self->{mrend} = $mrend;
}

=over 4

=item get_mcoords()

B<Description:> Retrieves the mlend, mrend for the cDNA coordinates.

B<Parameters:> none

B<Returns:> ($mlend, $mrend)

=back

=cut

sub get_mcoords () {
    my $self = shift;
    return ($self->{mlend}, $self->{mrend});
}


=over 4

=item get_per_id()

B<Description:> Retrieves the per_id for the alignment segment

B<Parameters:> none

B<Returns:> $per_id

=back

=cut

sub get_per_id {
    my $self = shift;
    return ($self->{per_id});
}



#private
sub set_orientation {
    my $self = shift;
    my $orientation = shift;
    $self->{orientation} = $orientation;
}

=over 4

=item get_orientation()

B<Description:> Retrieves the orientation for an alignment segment.

B<Parameters:> none

B<Returns:> [+|-]

=back

=cut

sub get_orientation {
    my $self = shift;
    return ($self->{orientation});
}

sub set_type {
    my $self = shift;
    my $type = shift;
    unless ($type =~ /first|last|internal|single/) {
        die "Incompatible segment type provided: $type\n";
    }
    $self->{type} = $type;
}


=over 4

=item get_type()

B<Description:>Retrieves the classification of the alignment segment 

B<Parameters:> none

B<Returns:> [first|last|internal|single]

=back

=cut

sub get_type {
    my $self = shift;
    return ($self->{type});
}

sub is_first {
    my $self = shift;
    return ($self->{type} eq "first") ;
}

sub is_internal {
    my $self = shift;
    return ($self->{type} eq "internal");
}

sub is_last {
    my $self = shift;
    return ($self->{type} eq "last");
}

sub is_single_segment {
    my $self = shift;
    return ($self->{type} eq "single");
}

sub set_left_splice_junction {
    my $self = shift;
    my $value = shift;
    $self->{has_left_splice_junction} = $value;
}

=over 4

=item has_left_splice_junction()

B<Description:> Provides result of a left splice junction test.

B<Parameters:> none

B<Returns:> [1|0]

1=true

0=false

=back

=cut


sub has_left_splice_junction {
    my $self = shift;
    return ($self->{has_left_splice_junction});
}

sub set_right_splice_junction {
    my $self = shift;
    my $value = shift;
    $self->{has_right_splice_junction} = $value;
}

=over 4

=item has_right_splice_junction()

B<Description:> Provides the result of a right splice junction test.

B<Parameters:> none

B<Returns:> [1|0]

=back

=cut


sub has_right_splice_junction {
    my $self = shift;
    return ($self->{has_right_splice_junction});
}

sub set_left_splice_site_chars () {
    my $self = shift;
    my $chars = shift;
    $self->{left_splice_site_chars} = $chars;
}


=over 4

=item get_left_splice_site_chars()

B<Description:> Retrieves the two characters representing the left splice site

B<Parameters:> none.

B<Returns:> $twochars

ie. Typically, this will return AG or AC depending on the spliced orientation.

=back

=cut

sub get_left_splice_site_chars () {
    my $self = shift;
    return ($self->{left_splice_site_chars});
}


sub set_right_splice_site_chars() {
    my $self = shift;
    my $chars = shift;
    $self->{right_splice_site_chars} = $chars;
}

=over 4

=item get_right_splice_chars()

B<Description:> Retrieves the two characters representing the right splice site

B<Parameters:> none.

B<Returns:> $two_chars

ie. typcially returns GT or CT depending on the spliced orientation.

=back

=cut


sub get_right_splice_site_chars() {
    my $self = shift;
    return ($self->{right_splice_site_chars});
}

sub toString() {
    my $self = shift;
    return( "segment\* orient: " . $self->{orientation} . " coords: " . $self->{lend} . "-" . $self->{rend} . " type:  " . $self->{type}
            . " lsplice: " . $self->has_left_splice_junction() . " rsplice: " . $self->has_right_splice_junction() . "\n");
    
}

sub toToken() {
    my $segment = shift;
    my $token = "";
    my $type = $segment->get_type();
    my ($lend, $rend) = $segment->get_coords();
    my ($mlend, $mrend) = $segment->get_mcoords();
    if ($type =~ /internal|last/) {
        # check splice site
        my $left_splice = $segment->get_left_splice_site_chars();
        if ($segment->has_left_splice_junction()) {
            $left_splice = uc $left_splice;
        } else {
            $left_splice = lc $left_splice;
        }
        $token .= $left_splice . "<";
    }
    $token .= $lend;
    if ($mlend) {
        $token .= "($mlend)";
    }
    $token .= "-$rend";
    if ($rend) {
        $token .= "($mrend)";
    }
    if ($type =~ /internal|first/) {
        # check splice site
        my $right_splice = $segment->get_right_splice_site_chars();
        if ($segment->has_right_splice_junction()) {
            $right_splice = uc $right_splice;
        } else {
            $right_splice = lc $right_splice;
        }
        $token .= ">$right_splice";
    }
    return ($token);
}



=over 4

=item clone()

B<Description:> Clones an Alignment_segment object into a new Alignment_segment object with the same attribute values.

B<Parameters:> none

B<Returns:> new CDNA::Alignment_segment

=back

=cut


sub clone {
    my $self = shift;
    my $packagename = ref $self;
    my $clone = {};
    foreach my $key (keys %$self) {
        $clone->{$key} = $self->{$key};
    }
    bless ($clone, $packagename);
    return ($clone);
}


=over 4

=item get_length()

B<Description:> Calculates the length of the alignment segment.

B<Parameters:> none

B<Returns:> int

=back

=cut



sub get_length {
    my $self = shift;
    my ($lend, $rend) = $self->get_coords();
    my $length = abs ($rend - $lend) + 1;
    return ($length);
}




1; #EOM
