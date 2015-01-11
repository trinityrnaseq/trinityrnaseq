package Simulate::Uniform_Read_Generator;

use strict;
use warnings;
use Carp;
use Data::Dumper;

sub new {
    my $packagename = shift;
    my $params_struct = shift;  
    
    # params should have format:
    # {  
    #  coordsets => [ [lend,rend], [lend,rend], ...],
    #  mean_fragment_length => int,
    #  fragment_length_stdev => int,
    #  read_length => int,
    # }
    
    unless (ref $params_struct eq "HASH") {
        confess "Error, params struct required";
    }
    
    
    my $coordsets_aref = $params_struct->{coordsets} or confess "Error, need coordsets parameter";
    my $mean_fragment_length = $params_struct->{mean_fragment_length} or confess "Error, need mean_fragment_length parameter";
    my $fragment_length_stdev = $params_struct->{fragment_length_stdev} or confess "Error, need fragment_length_stdev parameter";
    my $read_length = $params_struct->{read_length} or confess "Error, need read_length parameter";
           
    my $self = { _coordsets => undef,  # note, this is reconfigured below, not the same as the input parameter
                 mean_fragment_length => $mean_fragment_length,
                 fragment_length_stdev => $fragment_length_stdev,
                 read_length => $read_length,
             };
    
    bless($self, $packagename);
    
    $self->_init_coordsets($coordsets_aref);
    
    return($self);
}



####
sub _init_coordsets {
    my $self = shift;
    my ($coordsets_aref) = @_;

    my @ordered_coordsets = sort {$a->[0]<=>$b->[0]} @$coordsets_aref;
    
    my @coord_structs;

    my $prev_cdna_coord = 0;
    foreach my $coordset (@ordered_coordsets) {
        
        my ($lend, $rend) = @$coordset;

        my $cdna_lend = $prev_cdna_coord + 1;  ## cdna-relative  coordinates
        my $cdna_rend = $cdna_lend + ($rend - $lend);
        
        my $struct = { lend => $lend,
                       rend => $rend,
                       
                       cdna_lend => $cdna_lend,
                       cdna_rend => $cdna_rend,
        };

        push (@coord_structs, $struct);
    
        $prev_cdna_coord = $cdna_rend;
        
    }
    
    $self->{_coordsets} = \@coord_structs;

    return;
}


####
sub simulate_paired_reads_random_pos {
    my $self = shift;
    my ($num_reads) = @_;

    unless (defined($num_reads) && $num_reads =~ /\d/) {
        confess "Error, need num reads to simulate";
    }
    
    my $frag_length = $self->{mean_fragment_length};
    
    my $max_cdna_rend = $self->_get_max_cdna_rend();
    
    if ($frag_length > $max_cdna_rend) {
        $frag_length = $max_cdna_rend;
    }

        
    my @reads; # to contain (  [ left_read_aref, right_read_aref], ... )
    
    for (1..$num_reads) {
        
        my $read_start_pos = int( rand($max_cdna_rend - $frag_length + 1)) + 1;  # uniform selection of potential start sites.
        
        # simulate a fragment according to fragment length
        my @fragment = $self->_simulate_fragment($read_start_pos);
        
        # sample from each end to generate reads
        my @left_read_segs = $self->_get_left_fragment_read(@fragment);
        my @right_read_segs = $self->_get_right_fragment_read(@fragment);
        
        push (@reads, [  [@left_read_segs], [@right_read_segs] ] );
        
        ## simulate reads from fragment
        ## always do pairs and let the caller decide on which end or both to leverage.
        
    }
    
    
    return(@reads);
    
}


####
sub simulate_paired_reads_uniformly_across_seq {
    my $self = shift;
      
    my $frag_length = $self->{mean_fragment_length};
    
    my $max_cdna_rend = $self->_get_max_cdna_rend();
    
    if ($frag_length > $max_cdna_rend) {
        $frag_length = $max_cdna_rend;
    }
        
    my @reads; # to contain (  [ left_read_aref, right_read_aref], ... )
    
    for my $read_start_pos (1..($max_cdna_rend - $frag_length + 1)) {
                
        # simulate a fragment according to fragment length
        my @fragment = $self->_simulate_fragment($read_start_pos);
        
        # sample from each end to generate reads
        my @left_read_segs = $self->_get_left_fragment_read(@fragment);
        my @right_read_segs = $self->_get_right_fragment_read(@fragment);
        
        push (@reads, [  [@left_read_segs], [@right_read_segs] ] );
        
        ## simulate reads from fragment
        ## always do pairs and let the caller decide on which end or both to leverage.
        
    }
    
    
    return(@reads);
    
}



####
sub _simulate_fragment {
    my $self = shift;
    my ($read_start_pos) = @_;
    
    my $frag_length = $self->{mean_fragment_length};
    
    my $max_cdna_rend = $self->_get_max_cdna_rend();
    
    if ($frag_length > $max_cdna_rend) {
        $frag_length = $max_cdna_rend;
    }
            
    my @read_coords;
    
    my $structs_aref = $self->{_coordsets};
    
    my $read_cdna_lend_pos = $read_start_pos; #$self->_get_cdna_coord_via_genome_coord($read_start_pos);
    my $len_remaining = $frag_length;
    
    foreach my $struct (@$structs_aref) {
        
        my $exon_genome_lend = $struct->{lend};
        my $exon_genome_rend = $struct->{rend};

        my $exon_cdna_lend = $struct->{cdna_lend};
        my $exon_cdna_rend = $struct->{cdna_rend};

                
        if ($exon_cdna_lend <= $read_cdna_lend_pos && $read_cdna_lend_pos <= $exon_cdna_rend) {
            
            ## convert read lend position to genome coordinate position.
            my $delta = $read_cdna_lend_pos - $exon_cdna_lend;
            my $read_genome_lend = $exon_genome_lend + $delta;
            
            my $read_exon_len = $exon_cdna_rend - $read_cdna_lend_pos + 1;
            
            if ($read_exon_len >= $len_remaining) {
                my $read_genome_rend = $read_genome_lend + $len_remaining -1;
                push (@read_coords, [$read_genome_lend, $read_genome_rend]);
                last;
            }
            else {
                my $read_genome_rend = $exon_genome_rend;
                my $len_added = $read_genome_rend - $read_genome_lend + 1;
                push (@read_coords, [$read_genome_lend, $read_genome_rend]);
                $len_remaining -= $len_added;
                $read_cdna_lend_pos += $len_added;
            }
        }
        
        
    }
    
    return(@read_coords);
    
}




####
sub _get_max_cdna_rend {
    my $self = shift;

    my $structs_aref = $self->{_coordsets};
    
    my $max_cdna_rend = $structs_aref->[$#$structs_aref]->{cdna_rend};

    return($max_cdna_rend);
}


####
sub _get_cdna_coord_via_genome_coord {
    my $self = shift;
    my ($cdna_coord) = @_;
    
    my $structs_aref = $self->{_coordsets};
    
    foreach my $struct (@$structs_aref) {
        my $exon_genome_lend = $struct->{lend};
        my $exon_genome_rend = $struct->{rend};

        my $exon_cdna_lend = $struct->{cdna_lend};
        my $exon_cdna_rend = $struct->{cdna_rend};

        if ($cdna_coord >= $exon_cdna_lend && $cdna_coord <= $exon_cdna_rend) {
            
            my $delta = $cdna_coord - $exon_cdna_lend;
            my $genome_coord = $exon_genome_lend += $delta;
            return($genome_coord);
        }
    }

    confess "Error, could not map coordinate $cdna_coord within coordsets: " . Dumper($structs_aref);
}



####
sub _get_left_fragment_read {
    my $self = shift;
    my @frag = @_;

    my $read_length = $self->{read_length};

    my $length_remaining = $read_length;

    my @read_coords;
    
    foreach my $segment (@frag) {
        my ($frag_lend, $frag_rend) = @$segment;
        my $seg_len = $frag_rend - $frag_lend + 1;
        
        if ($seg_len >= $length_remaining) {
            my $read_seg = [$frag_lend, $frag_lend + $length_remaining - 1];
            push (@read_coords, $read_seg);
            last;
        }
        else {
            push (@read_coords, [$frag_lend, $frag_rend]);
            $length_remaining -= $seg_len;
        }
    }

    return(@read_coords);
}


####
sub _get_right_fragment_read {
    my $self = shift;
    my @frag = @_;

    my $read_length = $self->{read_length};
    
    my $length_remaining = $read_length;

    my @read_coords;

    foreach my $segment (reverse @frag) {
        my ($frag_lend, $frag_rend) = @$segment;
        my $seg_len = $frag_rend - $frag_lend + 1;
        
        if ($seg_len >= $length_remaining) {
            my $read_seg = [$frag_rend - $length_remaining + 1, $frag_rend];
            push (@read_coords, $read_seg);
            last;
        }
        else {
            push (@read_coords, [$frag_lend, $frag_rend]);
            $length_remaining -= $seg_len;
        }
    }


    @read_coords = reverse @read_coords;

    return(@read_coords);
}




1; #EOM
