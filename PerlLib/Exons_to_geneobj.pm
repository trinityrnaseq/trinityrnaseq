package main;
our $SEE;

package Exons_to_geneobj;

use strict;
use Longest_orf;
use Gene_obj;

## No reason to instantiate.  Use methods, fully qualified.
## allow for partial ORFS (missing start or stop codon in longest ORF)

###################
## Public method ##
###################

sub create_gene_obj {
    my ($exons_href, $sequence_ref, $partial_info_href) = @_;
    unless (ref $sequence_ref) {
        die "Error, need reference to sequence as input parameter.\n";
    }
    unless (ref $partial_info_href) {
        $partial_info_href = {};
    }
    
    ## exons_ref should be end5's keyed to end3's for all exons.
    my ($gene_struct_mod, $cdna_seq)  = &get_cdna_seq ($exons_href, $sequence_ref);
    
    my $cdna_seq_length = length $cdna_seq;
    my $long_orf_obj = new Longest_orf();
    
    # establish long orf finding parameters.
    $long_orf_obj->forward_strand_only();
    if ($partial_info_href->{"5prime"}) {
        print "Exons_to_geneobj: Allowing 5' partials\n" if $SEE;
        $long_orf_obj->allow_5prime_partials();
    }
    $long_orf_obj->allow_3prime_partials(); ## Allow this by default.
    
    
    $long_orf_obj->get_longest_orf($cdna_seq);
    my ($end5, $end3) = $long_orf_obj->get_end5_end3(); 
    
    print "CDS: $end5, $end3\n" if $SEE;
    my $gene_obj = &create_gene ($gene_struct_mod, $end5, $end3);
    
    
    ## Make adjustments if partial
    my @exons = $gene_obj->get_exons();
    my $first_exon = $exons[0];
    my $last_exon = $exons[$#exons];
    
    my $need_refine_flag = 0;
    
    if ($partial_info_href->{"5prime"}) {
        $gene_obj->set_5prime_partial(1);
        my $cds_exon = $first_exon->get_CDS_obj();
        if (ref $cds_exon) {
            $cds_exon->{phase} = 0;
            my $diff_coords = abs ($cds_exon->{end5} - $first_exon->{end5});
            if ($diff_coords > 0 && $diff_coords < 3) {
                $cds_exon->{end5} = $first_exon->{end5};
                $need_refine_flag = 1;
                
                ## update phase of first exon
                if ($diff_coords == 1) {
                    $cds_exon->{phase} = 2;
                } elsif ($diff_coords == 2) {
                    $cds_exon->{phase} = 1;
                }
                
            }
        }
        
    }
    
    if ($partial_info_href->{"3prime"}) {
        $gene_obj->set_3prime_partial(1);
        my $cds_exon = $last_exon->get_CDS_obj();
        if (ref $cds_exon) {
            my $diff_coords = abs ($cds_exon->{end3} - $last_exon->{end3});
            if ($diff_coords > 0 && $diff_coords < 3) {
                $cds_exon->{end3} = $last_exon->{end3};
                $need_refine_flag = 1;
            }
        }
    }
    
    ## Set the phase attribute for each cds exon:
    my @cds_exons;
    foreach my $exon ($gene_obj->get_exons()) {
        if (my $cds = $exon->get_CDS_obj()) {
            if (ref $cds) {
                push (@cds_exons, $cds);
            }
        }
    }
    
    if (@cds_exons) {
        my $cds_length = 0;
        my $first_cds = shift @cds_exons;
        my $phase = $first_cds->{phase};
        unless ($phase) {
            ## didn't set it for non 5' partials yet.
            $phase = $first_cds->{phase} = 0;
        }
        $cds_length -= $phase;
        $cds_length += $first_cds->length();
        foreach my $following_cds (@cds_exons) {
            $following_cds->{phase} = $cds_length % 3;
            $cds_length += $following_cds->length();
        }
    }
    
    
    if ($need_refine_flag) {
        $gene_obj->refine_gene_object();
    }
    
    return ($gene_obj);
}


########################
## Private methods #####
########################

####
sub get_cdna_seq {
    my ($gene_struct, $assembly_seq_ref) = @_;
    my (@end5s) = sort {$a<=>$b} keys %$gene_struct;
    my $strand = "?";
    foreach my $end5 (@end5s) {
	my $end3 = $gene_struct->{$end5};
	if ($end5 == $end3) { next;}
	$strand = ($end5 < $end3) ? '+':'-';
	last;
    }
    if ($strand eq "?") {
	print Dumper ($gene_struct);
	die "ERROR: I can't determine what orientation the cDNA is in!\n";
    }
    print NOTES "strand: $strand\n";
    my $cdna_seq;
    my $gene_struct_mod = {strand=>$strand,
		       exons=>[]}; #ordered lend->rend coordinate listing.
    foreach my $end5 (@end5s) {
	#print $end5;
	my $end3 = $gene_struct->{$end5};
	my ($coord1, $coord2) = sort {$a<=>$b} ($end5, $end3);
	my $exon_seq = substr ($$assembly_seq_ref, $coord1 - 1, ($coord2 - $coord1 + 1));
	$cdna_seq .= $exon_seq;
	push (@{$gene_struct_mod->{exons}}, [$coord1, $coord2]);
    }
    if ($strand eq '-') {
	$cdna_seq = reverse_complement ($cdna_seq);
    }
    return ($gene_struct_mod, $cdna_seq);
}


####
sub create_gene {
    my ($gene_struct_mod, $cds_pointer_lend, $cds_pointer_rend) = @_;
    my $strand = $gene_struct_mod->{strand};
    my @exons = @{$gene_struct_mod->{exons}};
    if ($strand eq '-') {
	@exons = reverse (@exons);
    }
    my $mRNA_pointer_lend = 1;
    my $mRNA_pointer_rend = 0;
    my $gene_obj = new Gene_obj();
    foreach my $coordset_ref (@exons) {
	my ($coord1, $coord2) = @$coordset_ref;
	my ($end5, $end3) = ($strand eq '+') ? ($coord1, $coord2) : ($coord2, $coord1);
	my $exon_obj = new mRNA_exon_obj($end5, $end3);
	my $exon_length = ($coord2 - $coord1 + 1);
	$mRNA_pointer_rend = $mRNA_pointer_lend + $exon_length - 1;
	## see if cds is within current cDNA range.
	if ( $cds_pointer_rend >= $mRNA_pointer_lend && $cds_pointer_lend <= $mRNA_pointer_rend) { #overlap
	    my $diff = $cds_pointer_lend - $mRNA_pointer_lend;
	    my $delta_lend = ($diff >0) ? $diff : 0;
	    $diff = $mRNA_pointer_rend - $cds_pointer_rend;
	    my $delta_rend = ($diff > 0) ? $diff : 0;
	    if ($strand eq '+') {
		$exon_obj->add_CDS_exon_obj($end5 + $delta_lend, $end3 - $delta_rend);
	    } else {
		$exon_obj->add_CDS_exon_obj($end5 - $delta_lend, $end3 + $delta_rend);
	    }
	}
	$gene_obj->add_mRNA_exon_obj($exon_obj);
	$mRNA_pointer_lend = $mRNA_pointer_rend + 1;
    }
    $gene_obj->refine_gene_object();
    #$gene_obj->{strand} = $strand;
    print $gene_obj->toString() if $SEE;
    return ($gene_obj);
}


sub reverse_complement { 
     my($s) = @_;
     my ($rc);
     $rc = reverse ($s);
     $rc =~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/;
     return($rc);
 }


	
1; #EOM
