#!/usr/local/bin/perl

package main;
our $SEE;


package CDNA::Alternative_splice_comparer;
use Gene_obj;
use strict;
use Data::Dumper;
use Carp;
use CDNA::PASA_alignment_assembler;

sub new {
    my $packagename = shift;
    my $self = {
        unspliced_introns => 0,
        conventional_alt_splice => 0,
        start_or_end_within_intron => 0,
        exon_skipping => 0,
        alternate_exons => 0
        };
    bless ($self, $packagename);
    return ($self);
}



####
sub compare_isoforms_via_alignmentObjs {
    my $self = shift;
    my ($align1, $align2) = @_;
    my $gene_1 = $align1->get_gene_obj_via_alignment();
    my $gene_2 = $align2->get_gene_obj_via_alignment();
    
    return ($self->compare_isoforms_via_geneObjs($gene_1, $gene_2));
}


####
sub compare_isoforms_via_geneObjs {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    ## Looking for:
    #    -unspliced introns
    #    -conventional alt-splice isoforms
    #    -transcriptional start or polyadenylation site within intron
    #    -exon skipping
    #    -alternate exons

    ## Look for Unspliced Introns
    my $struct  = { unspliced_introns => 0,
                    conventional_alt_splice => 0,
                    start_or_end_within_intron => 0,
                    exon_skipping => 0,
                    alternate_exons=> 0 };
    
    
    my @unspliced_introns = ($self->find_unspliced_introns($gene1, $gene2), $self->find_unspliced_introns($gene2, $gene1));
    if (@unspliced_introns) {
        print "*** Unspliced introns \n";
        $struct->{unspliced_introns} = 1;
        $self->{unspliced_introns} = \@unspliced_introns;
    }
    
    
    ## Look for the conventional alt splice isoforms (diff donors/acceptors for introns).
    my (%alternate_acceptors_n_donors) = $self->find_conventional_alt_splice_isoforms($gene1, $gene2);
    if (%alternate_acceptors_n_donors) {
        print "*** Conventional Alt splice (donor and/or acceptor)\n";
        #print Dumper (\%alternate_acceptors_n_donors);
        
        $self->{conventional_alt_splice} = \%alternate_acceptors_n_donors;
        if (@{$alternate_acceptors_n_donors{acceptors}}) {
            $struct->{conventional_alt_acceptor} = 1;
        }
        if (@{$alternate_acceptors_n_donors{donors}}) {
            $struct->{conventional_alt_donor} = 1;
        }
    }
    
    
    my %intron_starts_or_ends = ($self->find_starts_and_ends_within_introns ($gene1, $gene2), $self->find_starts_and_ends_within_introns ($gene2, $gene1));
    if (%intron_starts_or_ends) {
        print "*** Transcriptional start or polyadenylation site within an intron.\n";
        $struct->{start_or_end_within_intron} = 1;
        my @starts_or_ends = values %intron_starts_or_ends;
        $self->{start_or_end_within_intron} = \@starts_or_ends;
    }
    
    ## Look for exon skipping events.
    my @exon_skips = ($self->find_exon_skipping_events($gene1, $gene2), $self->find_exon_skipping_events($gene2, $gene1));
    if (@exon_skips) {
        print "*** Exon skipping event detected.\n";
        $struct->{exon_skipping} = 1;
        $self->{exon_skipping} = \@exon_skips;
    }
    
    ## Look or Alternate exons
    my @alternate_exons = ($self->find_alternate_exons($gene1, $gene2), $self->find_alternate_exons($gene2, $gene1));
    if (@alternate_exons) {
        print "*** Found alternate exons\n";
        $struct->{alternate_exons} = 1;
        $self->{alternate_exons} = \@alternate_exons;
    }
    
    return ($struct);
    
}


# private
sub enumerate_exons_of_gene {
    my $gene_obj = shift;
    # put everything in forward coordinate axis:
    my %exon_coords;
    my @exons = $gene_obj->get_exons();
    foreach my $exon (@exons) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        $exon_coords{$lend} = $rend;
    }
    return (%exon_coords);
}


#private
####
sub enumerate_introns_of_gene {
    my $gene_obj = shift;
    ## Put everything in forward strand coordinate axis.
    my %introns;
    my @exons = sort {$a->{end5}<=>$b->{end5}} $gene_obj->get_exons();
    for (my $i = 0; $i < $#exons; $i++) {
        my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $exons[$i]->get_coords();
        my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $exons[$i+1]->get_coords();
        my ($intron_end5, $intron_end3) = ($exon1_rend + 1, $exon2_lend - 1);
        $introns{$intron_end5} = $intron_end3;
    }
    return (%introns);
}


=over 4

=item find_unspliced_introns()

B<Description:> Find unspliced introns in gene_1 when compared to gene_2

B<Parameters:> $gene1, $gene2

B<Returns:> @unspliced_introns

@unspliced_introns is a list of coordinate pairs representing the unspliced introns found in gene 1 when compared to gene2

@unspliced_introns = ([end5,end3], ...)

=back

=cut

####
sub find_unspliced_introns {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    ## Look for unspliced intron found in gene1 when compared to gene2
    my %gene1_exon_coords = &enumerate_exons_of_gene ($gene1);
    my %gene2_intron_coords = &enumerate_introns_of_gene($gene2);
    
    my @unspliced_introns;
    foreach my $intron_lend (keys %gene2_intron_coords) {
        my $intron_rend = $gene2_intron_coords{$intron_lend};
        
        foreach my $exon_lend (keys %gene1_exon_coords) {
            my $exon_rend = $gene1_exon_coords{$exon_lend};
            
            if ($intron_lend > $exon_lend && $intron_rend < $exon_rend) { #unspliced intron found
                push (@unspliced_introns, [$intron_lend, $intron_rend]);
            }
        }
    }
    return (@unspliced_introns);
}



=over 4

=item find_conventional_alt_splice_isoforms()

B<Description:> Looks for different donor and acceptor sites within overlapping introns of genes

B<Parameters:> $gene1, $gene2

B<Returns:> %alt_donors_and_acceptors

with structure:

%alt_donors_and_acceptors = ( acceptors => 
                              [    
                                   { gene1 => acceptor_coord, gene2 => acceptor_coord }, ...
                                   
                                   
                                   
                                   ],
                              
                              donors => [
                                         
                                         { gene1 => donor_coord, gene2 => donor_coord }, ...
                                         
                                         
                                         ]
                              
                              );

     Coordinates stored are the actual exon boundary coordinates (first or last bp of each exon)

    


=back

=cut


####
sub find_conventional_alt_splice_isoforms {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    print "## Looking for conventional alt splice isoforms (diff donors, acceptors)\n" if $SEE;
    my %exons_1_hash = &enumerate_exons_of_gene ($gene1);
    my %exons_2_hash = &enumerate_exons_of_gene ($gene2);
    
    my $orientation = $gene1->get_orientation();
    if ($orientation ne $gene2->get_orientation()) {
        die "Error, inconsistent orientations between genes: " . $gene1->toString() . $gene2->toString();
    }
    
    # algorithm
    #   -find one-to-one mappings between exons, and locate differences at acceptors and donor sites.
    
    my %alternate_acceptors_n_donors = ( acceptors => [],
                                         donors => [] 
                                         ); # holds coordinates for all gene1 diff boundaries.
    
    my $found_diff_flag = 0;
    
    
    ## compare exons of gene1 to gene2
    my @exons_gene_1_list;
    my @exons_gene_2_list;
    # build data structure:
    foreach my $data_pair ( [\%exons_1_hash, \@exons_gene_1_list], 
                            [\%exons_2_hash, \@exons_gene_2_list] ) {
        
        my ($exons_href, $exons_aref) = @$data_pair;
        
        foreach my $lend (keys %$exons_href) {
            my $rend = $exons_href->{$lend};
            
            push (@$exons_aref, { lend => $lend,
                                  rend => $rend,
                                  match_indices => [] } );
        }
    }

    @exons_gene_1_list = sort {$a->{lend}<=>$b->{lend}} @exons_gene_1_list;
    @exons_gene_2_list = sort {$a->{lend}<=>$b->{lend}} @exons_gene_2_list;
    

    # all-vs-all comparison:
    for (my $i = 0; $i <= $#exons_gene_1_list; $i++) {

        my $i_ele_ref = $exons_gene_1_list[$i];
        my ($i_lend, $i_rend, $i_match_indices_aref) = ($i_ele_ref->{lend},
                                                        $i_ele_ref->{rend},
                                                        $i_ele_ref->{match_indices} );
        

        for (my $j = 0; $j <= $#exons_gene_2_list; $j++) {

            my $j_ele_ref = $exons_gene_2_list[$j];
            my ($j_lend, $j_rend, $j_match_indices_aref) = ($j_ele_ref->{lend},
                                                            $j_ele_ref->{rend},
                                                            $j_ele_ref->{match_indices});
            
            if ($i_lend < $j_rend && $i_rend > $j_lend) { #overlap
                push (@$i_match_indices_aref, $j);
                push (@$j_match_indices_aref, $i);
            }
        }
    }

    ## find donors and acceptors:
    ## check gene_1's exons for 1-1 mappings and end differences at splice junctions
    
    for (my $i = 0; $i <= $#exons_gene_1_list; $i++) {
        
        my $i_ele_ref = $exons_gene_1_list[$i];
        my ($i_lend, $i_rend, $i_match_indices_aref) = ($i_ele_ref->{lend},
                                                        $i_ele_ref->{rend},
                                                        $i_ele_ref->{match_indices} );
        
        if (scalar (@$i_match_indices_aref) == 1) {
            ## found some mapping
            my $j_index = $i_match_indices_aref->[0];
            my $j_ele_ref = $exons_gene_2_list[$j_index];
            my ($j_lend, $j_rend, $j_match_indices_aref) = ($j_ele_ref->{lend},
                                                            $j_ele_ref->{rend},
                                                            $j_ele_ref->{match_indices});
            

            if (scalar (@$j_match_indices_aref) != 1) {
                next; ## this is a 1-many mapping, want only 1-1 mappings
            }

            # make sure j's i is i
            if ($j_match_indices_aref->[0] != $i) {
                ## bad, this should never happen!
                confess "Error, found exon 1-1 mapping of $i to $j_index, but j maps to @$j_match_indices_aref ";
            }
            

            ## check left boundary:
            if ($i_lend != $j_lend ## diff coordinate
                && $i != 0  # at a splice junction
                && $j_index != 0  # at a splice junction
                ) {
                
                ## found splice difference at left junction:
                $found_diff_flag = 1;
                
                my $splice_ref = ($orientation eq '+') 
                    ? $alternate_acceptors_n_donors{acceptors}
                    : $alternate_acceptors_n_donors{donors};
                
                push (@$splice_ref, { gene1 => $i_lend,
                                      gene2 => $j_lend } );
                
            }
            

            ## check right boundary:
            if ($i_rend != $j_rend ## diff coordinate
                && $i != $#exons_gene_1_list # at splice junction
                && $j_index != $#exons_gene_2_list # at splice junction
                ) {
                
                $found_diff_flag = 1;
                
                my $splice_ref = ($orientation eq '+') 
                    ? $alternate_acceptors_n_donors{donors}
                    : $alternate_acceptors_n_donors{acceptors};
                
                push (@$splice_ref, { gene1 => $i_rend,
                                      gene2 => $j_rend } );
            }
        }
    }
    
    if ($found_diff_flag) {
        return (%alternate_acceptors_n_donors);
    } else {
        return ();
    }
    
}




=over 4

=item find_exon_skipping_events()

B<Description:> Finds an exon of gene_1 which reside within an intron of gene_2

B<Parameters:> gene1, gene2

B<Returns:> @skipped_exons


 notice this is a list of lists
 each list is a set of adjacent skipped exons, joined so that they correspond to a single event.
so what we are really getting here is a list of events of skipped exons where each event may contain one or more skipped exons.


@skipped_exons = ( 

                   [
                    [exon_lend,exon_rend], ...
                    
                    ],

                   [
                    [exon_lend, exon_rend], ...

                    ]

                   
                    ) 

=back

=cut


####
sub find_exon_skipping_events {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    
    # Algorithm:
    #   -find an internal exon of gene 1 that resides within an intron of gene 2.  Flanking exons must be anchored to the other isoform
    my %gene1_exons = &enumerate_exons_of_gene($gene1);
    my %gene2_introns = &enumerate_introns_of_gene($gene2);
    
    my @potential_skipped_exons;
    foreach my $exon1_lend (keys %gene1_exons) {
        my $exon1_rend = $gene1_exons{$exon1_lend};
        
        ## See if within intron of second gene
        foreach my $intron2_lend (keys %gene2_introns) {
            my $intron2_rend = $gene2_introns{$intron2_lend};
            
            if ($exon1_lend > $intron2_lend && $exon1_rend < $intron2_rend) { #exon incapsulated in intron
                push (@potential_skipped_exons, [$exon1_lend, $exon1_rend]);
            }
        }
    }
    
    ## Verify flanking exons are anchorable:
    my @skipped_exons;
    if (@potential_skipped_exons) {
        foreach my $potential_skipped_exon (@potential_skipped_exons) {
            my ($exon_lend, $exon_rend) = @$potential_skipped_exon;
            
            ## Try to anchor left exon
            my $anchor_left_exon = 0;
            foreach my $exon1 ($gene1->get_exons()) {
                my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $exon1->get_coords();
                unless ($exon1_rend < $exon_lend) { next;}
                foreach my $exon2 ($gene2->get_exons()) {
                    my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $exon2->get_coords();
                    unless ($exon2_rend < $exon_lend) { next;}
                    
                    if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) { #anchorable
                        $anchor_left_exon = 1;
                        last;
                    }
                }
                if ($anchor_left_exon) { last;}
            }
            unless ($anchor_left_exon) { next;}
            
            ## Try to anchor the right exon
            my $anchor_right_exon = 0;
            foreach my $exon1 ($gene1->get_exons()) {
                my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $exon1->get_coords();
                unless ($exon1_lend > $exon_rend) { next;}
                foreach my $exon2 ($gene2->get_exons()) {
                    my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $exon2->get_coords();
                    unless ($exon2_lend > $exon_rend) { next;}
                    
                    if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) { #anchorable
                        $anchor_right_exon = 1;
                        last;
                    }
                }
                if ($anchor_right_exon) { last;}
            }
            if ($anchor_right_exon && $anchor_left_exon) {
                push (@skipped_exons, $potential_skipped_exon);
            }
        }
        
    }


    ## group into lists of adjacent exons
    
    my @ret_skipped_exons;
    if (@skipped_exons) {
        @skipped_exons = sort {$a->[0]<=>$b->[0]} @skipped_exons;
        
        my %coord_to_order;
        ## map each exon to an integer
        my $order = 0;
        foreach my $exon (sort {$a->{end5}<=>$b->{end5}} $gene1->get_exons()) {
            my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
            $order++;
            $coord_to_order{$lend} = $order;
        }
        
        my $first_skipped_exon = shift @skipped_exons;
        @ret_skipped_exons = ([$first_skipped_exon]);
        
        while (@skipped_exons) {
            my $last_event = $ret_skipped_exons[$#ret_skipped_exons];
            my $last_skipped_exon = $last_event->[ $#{$last_event} ];
            
            my $last_lend = $last_skipped_exon->[0];
            
            my $curr_skipped_exon = shift @skipped_exons;
            my $curr_lend = $curr_skipped_exon->[0];
            
            if ($coord_to_order{$curr_lend} - $coord_to_order{$last_lend} == 1) {
                ## adjacent, so group them
                push (@$last_event, $curr_skipped_exon);
            }
            else {
                ## not adjacent
                # start new event
                push (@ret_skipped_exons, [$curr_skipped_exon]);
            }
        }
    }
    
    return (@ret_skipped_exons);
    
}







=over 4

=item find_alternate_exons()

B<Description:> Finds terminal exons in gene1 that are different and non-overlapping, and adjacent to overlapping exons.

B<Parameters:> $gene1, $gene2

B<Returns:> @range_of_coords_containing_alternate_exons

    @ret = ( {  type => lend|rend,
                coords => [region_lend,region_rend],
                num_exons => intval
                }
             
             , ...
             
             )




=back

=cut

sub find_alternate_exons {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    my @alternate_exon_regions; # store coords of alternate exons
    # Algorithm:
    #    -Looking at terminal exons, should have non-overlapping exons prior to the first overlapping exon
    
    ## Look from front to back:
    my @gene1_exons = sort {$a->{end5}<=>$b->{end5}} $gene1->get_exons();
    my @gene2_exons = sort {$a->{end5}<=>$b->{end5}} $gene2->get_exons();
    
    my @alternate_exons_front;
    for (my $i = 0; $i <= $#gene1_exons; $i++) {
        my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $gene1_exons[$i]->get_coords();
        my $overlapping_j = undef();
        for (my $j = 0; $j <= $#gene2_exons; $j++) {
            my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $gene2_exons[$j]->get_coords();
            
            ## check for overlap
            if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) {
                $overlapping_j = $j;
                last;
            }
        }
        if (defined($overlapping_j)) {
            ## See if i and j are not first:
            if ($i != 0 && $overlapping_j != 0) {
                for (my $x=0; $x < $i; $x++) {
                    my $exon = $gene1_exons[$x];
                    my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
                    push (@alternate_exons_front, [$lend,$rend]);
                }
            }
            last;
        }
    }
    
    ## Look from back to front:

    my @alternate_exons_back;
    for (my $i = $#gene1_exons; $i >= 0; $i--) {
        my ($exon1_lend, $exon1_rend) = sort {$a<=>$b} $gene1_exons[$i]->get_coords();
        my $overlapping_j = undef();
        for (my $j = $#gene2_exons; $j >= 0; $j--) {
            my ($exon2_lend, $exon2_rend) = sort {$a<=>$b} $gene2_exons[$j]->get_coords();
            
            ## check for overlap
            if ($exon1_lend < $exon2_rend && $exon1_rend > $exon2_lend) {
                $overlapping_j = $j;
                last;
            }
        }
        if (defined($overlapping_j)) {
            ## See if i and j are not last:
            if ($i != $#gene1_exons && $overlapping_j != $#gene2_exons) {
                for (my $x=$#gene1_exons; $x > $i; $x--) {
                    my $exon = $gene1_exons[$x];
                    my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
                    push (@alternate_exons_back, [$lend,$rend]);
                }
            }
            last;
        }
    }
    
    if (@alternate_exons_front) {
        
        my @front_coords;
        my $num_alternate_exons_front = scalar (@alternate_exons_front);
        foreach my $coordpair (@alternate_exons_front) {
            push (@front_coords, @$coordpair);
        }
        @front_coords = sort {$a<=>$b} @front_coords;
        my $region_lend = shift @front_coords;
        my $region_rend = pop @front_coords;
        push (@alternate_exon_regions, { type => 'lend',
                                         coords => [$region_lend, $region_rend],
                                         num_exons => $num_alternate_exons_front,
                                         }
              );
    }

    if (@alternate_exons_back) {
        my @back_coords;
        my $num_alternate_exons_back = scalar (@alternate_exons_back);
        foreach my $coordpair (@alternate_exons_back) {
            push (@back_coords, @$coordpair);
        }
        @back_coords = sort {$a<=>$b} @back_coords;
        my $region_lend = shift @back_coords;
        my $region_rend = pop @back_coords;
        push (@alternate_exon_regions, { type => 'rend',
                                         coords => [$region_lend, $region_rend],
                                         num_exons => $num_alternate_exons_back
                                             }
              );
    }
    

    return (@alternate_exon_regions);
}


=over 4

=item find_starts_and_ends_within_introns()

B<Description:> The first and last exons of gene_1 are compared to the introns of gene_2. 

B<Parameters:> $gene1, $gene2

B<Returns:> @coords

@coords contains the coordinates of either the very end5 or very end3 of terminal exons which fall into introns of gene_2

=back

=cut


####
sub find_starts_and_ends_within_introns {
    my $self = shift;
    my ($gene1, $gene2) = @_;

    my $fuzzlength = $CDNA::PASA_alignment_assembler::FUZZLENGTH;
    
    
    my %starts_and_ends;

    my $orientation = $gene1->get_orientation();

    # Algorithm:
    #   -first and last exon of gene1 is compared to introns of gene2
    my @gene1_exons = $gene1->get_exons();
    my %gene2_introns = &enumerate_introns_of_gene($gene2);
    my %gene2_exons = &enumerate_exons_of_gene($gene2);
    my $first_exon = $gene1_exons[0];
    my ($end5, $end3) = $first_exon->get_coords();

    foreach my $intron_lend (keys %gene2_introns) {
        my $intron_rend = $gene2_introns{$intron_lend};
        if ($end5 >= $intron_lend && $end5 <= $intron_rend) { #endpoint encapsulated by intron.
    
            ## make sure it's not fuzz:
            if ($orientation eq "+") {
                if ( abs ($end5-$intron_rend) + 1 <= $fuzzlength) {
                    next;
                }
            }
            else { # minus strand
                if (abs ($end5-$intron_lend)+1 <= $fuzzlength) {
                    next;
                }
            }
            
            ## make sure exon overlaps another exon
            foreach my $exon_lend (keys %gene2_exons) {
                my $exon_rend = $gene2_exons{$exon_lend};
                my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
                if ($rend > $exon_lend && $lend < $exon_rend) { #overlap
                    $starts_and_ends{start} = $end5; # store start
                    last;
                }
            }
            last;
        }
    }
    
    ## Now try last exon
    my $last_exon = $gene1_exons[$#gene1_exons];
    my ($end5, $end3) = $last_exon->get_coords();

    foreach my $intron_lend (keys %gene2_introns) {
        my $intron_rend = $gene2_introns{$intron_lend};
        if ($end3 >= $intron_lend && $end3 <= $intron_rend) { #endpoint encapsulated by intron.
            
            ## make sure not fuzz:
            if ($orientation eq "+") {
                if (abs ($end3 - $intron_lend) + 1 <= $fuzzlength) {
                    next;
                }
            }
            else { #minus strand
                if (abs ($end3 - $intron_rend) + 1 <= $fuzzlength) {
                    next;
                }
            }
            
            ## Make sure exon overlaps another exon
            foreach my $exon_lend (keys %gene2_exons) {
                my $exon_rend = $gene2_exons{$exon_lend};
                my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
                if ($rend > $exon_lend && $lend < $exon_rend) {
                    $starts_and_ends{end} = $end3; #store end
                    last;
                }
            }
            last;
        }
    }
    return (%starts_and_ends);
}



=over 4

=item compare_exons()

B<Description:> Compares all CDS exons between genes 1 and 2, returns number of identical CDS exons and total number of CDS exons between the two genes.

B<Parameters:> $gene1, $gene2

B<Returns:> ($num_identical_CDS_exons, $num_total_CDS_exons)


=back

=cut



####
sub compare_exons {
    my $self = shift;
    my ($gene1, $gene2) = @_;
    print "gene1_strand: $gene1->{strand}\n";
    my $clone_1 = $gene1->clone_gene();
    print "clone1_strand: " . $clone_1->{strand} . "\n";
    my $clone_2 = $gene2->clone_gene();
    $clone_1->trim_UTRs();
    print "clone1_strand, utrs trimmed: " . $clone_1->{strand} . "\n";
    $clone_2->trim_UTRs();
   
    my @exons_1 = $clone_1->get_exons();
    
    my @exons_2 = $clone_2->get_exons();
    
    my @identity_list = ();
    my @all_exons = sort {$a->{end5}<=>$b->{end5}} (@exons_1, @exons_2);
    
    for (my $i=0; $i <= $#all_exons-1; $i++) {
        
        my $curr_exon = $all_exons[$i];
        my $next_exon = $all_exons[$i+1];
        
        my ($curr_exon_end5, $curr_exon_end3) = $curr_exon->get_coords();
        
        my ($next_exon_end5, $next_exon_end3) = $next_exon->get_coords();
        
        if ($curr_exon_end5 == $next_exon_end5 && $curr_exon_end3 == $next_exon_end3) {
            $identity_list[$i] = 1;
            $identity_list[$i+1] = 1;
            $i++; #if A = B, then go onto comparing C to D, not B to C.
        }
    }
    
    my ($num_identical_exons, $total_num_exons) = (0,0);
    for (my $i=0; $i <= $#all_exons; $i++) {
        if ($identity_list[$i]) {
            $num_identical_exons++;
        }
        $total_num_exons++;
    }
    
    return ($num_identical_exons, $total_num_exons);
}

	


=over 4

=item start_or_stop_within_intron()

B<Description:> Compares the start codon and stop codon position of gene1 to the introns of gene2.

B<Parameters:> $gene1, $gene2

B<Returns:> ($start_within_intron, $stop_within_intron)

return values are 0|1 meaning true|false for each return parameter.

=back

=cut


sub start_or_stop_codon_within_intron {
    my $self = shift;
    my ($gene1, $gene2) = @_;
 

    ## Look for annotated start codon or stop codon within intron:
    my ($annotated_start_within_intron, $annotated_stop_within_intron) = (0,0);
    
    my ($start_codon, $stop_codon) = $gene1->get_model_span();
    my (@alignment_segments) = sort {$a->{end5}<=>$b->{end5}} $gene2->get_exons();
    if ($#alignment_segments > 0) { #multiple segments:
        for (my $i=1; $i <= $#alignment_segments; $i++) {
            my $prev_seg = $alignment_segments[$i-1];
            my ($prev_lend, $prev_rend) = sort {$a<=>$b} $prev_seg->get_coords();
            my $curr_seg = $alignment_segments[$i];
            my ($curr_lend, $curr_rend) = sort {$a<=>$b} $curr_seg->get_coords();
            my ($intron_lend, $intron_rend) = ($prev_rend+1, $curr_lend-1);
            if ($start_codon >= $intron_lend && $start_codon <= $intron_rend) {
                $annotated_start_within_intron = 1;
            }
            if ($stop_codon >= $intron_lend && $stop_codon <= $intron_rend) {
                $annotated_stop_within_intron = 1;
            }
        }
    }
    return ($annotated_start_within_intron, $annotated_stop_within_intron);
}






#####
## Static methods
####

=over 4

=item adjust_alternate_exon_region_coords()

B<Description:> Adjusts the coordinates provided by find_alternate_exons() so that they are extended to include the adjacent intron.  

B<Parameters:> gene_obj, region_lend, region_rend

B<Returns:> adjusted_region_lend, adjusted_region_rend


This is useful in cases where we try to see if the variation impacts the protein coding region when compared to an alternate gene.

region_lend and region_rend are the lend,rend values stored in table: splice_variation under type = 'alternate_exon'


=back

=cut

sub adjust_alternate_exon_region_coords {
    
    my ($gene_obj, $region_lend, $region_rend) = @_;
    
    my @other_exon_coords;
    foreach my $exon ($gene_obj->get_exons()) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        unless ($lend >= $region_lend && $rend <= $region_rend) {
            # not encapsulated:
            push (@other_exon_coords, $lend, $rend);
        }
    }
    @other_exon_coords = sort {$a<=>$b} @other_exon_coords;
    
    #print "other exon coords: @other_exon_coords\n";
    
    my $other_lend = shift @other_exon_coords;
    my $other_rend = pop @other_exon_coords;
    
    ## make sure there isn't any overlap
    if ($other_lend <= $region_rend && $other_rend >= $region_lend) {
        #overlap BAD!
        confess ("Error, coordinate regions overlap ($other_lend, $other_rend) w/ region($region_lend, $region_rend)\n"
                 . Dumper (\@other_exon_coords));
    }
    
    ## adjust boundary to include intron region
    if ($other_rend < $region_lend) {
        $region_lend = $other_rend + 1;
    }
    elsif ($other_lend > $region_rend) {
        $region_rend = $other_lend - 1;
    }
    else {
        confess "Error, cannot figure out how to adjust the boundaries." . Dumper (\@other_exon_coords);
    }
    
    return ($region_lend, $region_rend);
}


=over 4

=item extend_coords_to_intron_bounds()

B<Description:> given the coordinates of a region of skipped exons, the coordinates are extended to the far bounds of the adjacent introns

B<Parameters:> (gene_obj, region_lend, region_rend)

B<Returns:> (adjusted_region_lend, adjusted_region_rend)

=back

=cut


sub extend_coords_to_intron_bounds {
    my ($gene_obj, $region_lend, $region_rend) = @_;
    
    my @left_coords;
    my @right_coords;
    
    foreach my $exon ($gene_obj->get_exons()) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        
        if ($lend > $region_rend) {
            push (@right_coords, $lend, $rend);
        }
        elsif ($rend < $region_lend) {
            push (@left_coords, $lend, $rend);
        }
        
    }

    @left_coords = sort {$a<=>$b} @left_coords;
    @right_coords = sort {$a<=>$b} @right_coords;

    unless (@left_coords && @right_coords) {
        confess "Error, missing either left or right coords: \n"
            . "left: @left_coords\n"
            . "right: @right_coords\n"
            . "region: $region_lend, $region_rend\n";
    }
    

    my $new_left_bound = pop @left_coords;
    $new_left_bound++;
    
    my $new_right_bound = shift @right_coords;
    $new_right_bound--;
    
    return ($new_left_bound, $new_right_bound);
}





1;
	
    
