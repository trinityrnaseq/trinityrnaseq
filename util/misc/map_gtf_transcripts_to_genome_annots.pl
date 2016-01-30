#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib/");
use Gene_obj;
use GFF3_utils;
use GTF_utils;
use SAM_reader;
use SAM_entry;
use Overlap_info;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Data::Dumper;
use Overlap_piler;

$ENV{LC_ALL} = 'C';


my $usage = <<_EOUSAGE_;

#########################################################################
#
#  *Required:

#  # Annot settings:
#
#  --annot_gff3 <string>     gff3 file name
#    OR
#  --annot_gtf <string>      gtf file name 
#
#  # Transcript alignment settings:
#  
#  --alignment_gtf <string>   alignments in gtf format
#
# 
# *Optional:
#
#  # Intergenic settings: (default, off)
#
#  --include_intergenic          create features out of the intergenic regions.
#  --min_intergenic_length <int> minimum size of an intergenic feature to be included.
#
#  # Reporting settings:
#
#  --best                        only report the mapping that has the highest % overlap with the transcript
#  --ignore_antisense            only consider sense mappings  
#  --ignore_strandedness         ignore all strand orientation, set to '?' as done for intergenic regions.
#
#  --no_require_compatibility    do not require compatible overlap (default: does. splice junctions must match up)
#
###################################################################################################


_EOUSAGE_

	;


my $annot_genes_gff3;
my $annot_genes_gtf;

my $alignment_gtf;

my $help_flag;

my $VERBOSE = 0;
my $FUZZY_OVERLAP = 0;

my $INCLUDE_INTERGENIC = 0;


my $MIN_INTERGENIC_LENGTH = 100;
my $MAX_MERGE_INDEL = 5;

my $BEST_ONLY = 0;
my $IGNORE_ANTISENSE = 0;
my $IGNORE_STRANDEDNESS = 0;

my $NO_REQUIRE_COMPATIBILITY = 0;

&GetOptions ( 'h' => \$help_flag,

              'annot_gff3=s' => \$annot_genes_gff3,
              'annot_gtf=s' => \$annot_genes_gtf,
              
              'alignment_gtf=s' => \$alignment_gtf,
              
              
              'v' => \$VERBOSE,
              'FUZZY' => \$FUZZY_OVERLAP,

              'include_intergenic' => \$INCLUDE_INTERGENIC,
              'min_intergenic_length=i' => \$MIN_INTERGENIC_LENGTH,
              
              'best' => \$BEST_ONLY,
              'ignore_antisense' => \$IGNORE_ANTISENSE,
              'ignore_strandedness' => \$IGNORE_STRANDEDNESS,
              
              'no_require_compatibility' => \$NO_REQUIRE_COMPATIBILITY,
              );


if ($help_flag) { die $usage; }


unless ($annot_genes_gff3 || $annot_genes_gtf) { 
	die $usage;
}


if (@ARGV) {
    die "Options: @ARGV not understood";
}

main: {
	
    my $gene_obj_indexer_href = {};

	my $contig_to_gene_list_href;

    print STDERR "-parsing gene annotations file\n";
	if ($annot_genes_gff3) {
        
		## associate gene identifiers with contig id's.
		$contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($annot_genes_gff3, $gene_obj_indexer_href);
	}
	else {
		# GTF mode
		$contig_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($annot_genes_gtf, $gene_obj_indexer_href);
	}
    
        
    ## Populate Genomic Features

    print STDERR "-organizing annotated transcript features.\n";

	my %chr_to_features = &get_chr_to_features($contig_to_gene_list_href, $gene_obj_indexer_href);

    if ($INCLUDE_INTERGENIC) {
        print STDERR "-adding intergenic features\n";
        &add_intergenic_features(\%chr_to_features);
    }

	my %feature_lengths = &get_feature_lengths(\%chr_to_features);
    
    
    ## Populate Transcript Alignment Features
    my $trans_obj_indexer_href = {};
    my $contig_to_trans_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($alignment_gtf, $trans_obj_indexer_href);
    my %trans_align_features = &get_chr_to_features($contig_to_trans_list_href, $trans_obj_indexer_href);
    
    
    ## assign read mapping to features, including sense and antisense mappings
    print STDERR "\n\n-mapping transcripts to features\n";
    my $mapped_trans_file = &map_transcript_alignments_to_genome_features(\%chr_to_features, \%trans_align_features);
    
    
    
    print STDERR "\n\n-refining mappings\n\n";
    &refine_mappings_estimate_counts($mapped_trans_file); 
    

    exit(0);
    

}


####
sub get_chr_to_features {
    my ($contig_to_gene_list_href, $gene_obj_indexer_href) = @_;


    my %chr_to_features;
	
    my $total_transcripts = 0;

	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
		foreach my $gene_id (@gene_ids) {
			my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
			
			unless (ref $gene_obj_ref) {
				die "Error, no gene_obj for gene_id: $gene_id";
			}
			
			my $strand = $gene_obj_ref->get_orientation();

			
			my $scaffold = $asmbl_id;
            
			my $max_isoform_cdna_length = 0;
			
			foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                
                $total_transcripts++;
                
                my $isoform_id = join("::", $gene_id, $isoform->{Model_feat_name}); ## embed the gene identifier in with the transcript id, so we can tease it apart later.
				my $isoform_length = 0;
				
                my @coordset = &get_coordset_for_isoform($isoform);
                my ($lend, $rend) = ($coordset[0]->[0], $coordset[$#coordset]->[1]);
                
                my $length = &sum_coordset_segments(\@coordset);

                
                push (@{$chr_to_features{$scaffold}}, { 
                    acc => $isoform_id,
                    name => $isoform->{com_name},
                    lend => $lend,
                    rend => $rend,
                    coordset => [@coordset],
                    strand => $strand,
                    length => $length,
                      } );
                

            }
		}
	}
         
    return(%chr_to_features);
   
}



####
sub add_intergenic_features {
    my ($chr_to_features_href) = @_;

    my $inter_counter = 0;
    
    foreach my $feature_list_aref (values %$chr_to_features_href) {
        
        my @coords;
        foreach my $feature (@$feature_list_aref) {
            my ($lend, $rend) = ($feature->{lend}, $feature->{rend});
            push (@coords, [$lend,$rend]);
        }
        
        my @piles = &Overlap_piler::simple_coordsets_collapser(@coords);
        
        my $prev_rend = undef;
        foreach my $pile (@piles) {
            my ($curr_lend, $curr_rend) = @$pile;
            if ($prev_rend) {
                my $inter_lend = $prev_rend + 1;
                my $inter_rend = $curr_lend - 1;
                
                my $inter_len = $inter_rend - $inter_lend + 1;

                if ($inter_lend < $inter_rend && $inter_len >= $MIN_INTERGENIC_LENGTH) {
                    $inter_counter++;
                    
                    push (@{$feature_list_aref}, {
                        acc => "intergeneic_region.$inter_counter",
                        name => "intergenic region",
                        lend => $inter_lend,
                        rend => $inter_rend,
                        coordset => [ [$inter_lend, $inter_rend] ],
                        strand => '+',
                        length => $inter_len,
                        } );
                }
                
            }
            $prev_rend = $curr_rend;
        }

    }


    return;
}



####
sub map_transcript_alignments_to_genome_features {
    my ($chr_to_features_href, $trans_align_features_href) = @_;
    
	##
	## Map transcripts to genes
	##


	
    print STDERR "-examining transcript alignments, mapping to annotated features.\n";
    

    my $trans_mapping_file = "trans_align_mappings.$$.txt";
    open (my $ofh, ">$trans_mapping_file") or die $!;
    
    foreach my $scaffold (sort keys %$trans_align_features_href) {
        

        my @chr_features;
        if (exists $chr_to_features_href->{$scaffold}) { 
            @chr_features = sort {$a->{lend} <=> $b->{lend}} @{$chr_to_features_href->{$scaffold}};
            
        }
        
        my @trans_features = sort {$a->{rend}<=>$b->{rend}} @{$trans_align_features_href->{$scaffold}};
        
        foreach my $trans_feature (@trans_features) {
            my $position_lend = $trans_feature->{lend};
            my $position_rend = $trans_feature->{rend};
                        
            my $trans_feature_coordset = $trans_feature->{coordset};
            my $trans_feature_strand = $trans_feature->{strand};
            

            my @container;  ## holds current features.
            
            ## collect features that overlap 
            while (@chr_features && 
                   $chr_features[0]->{lend} <= $position_rend) {
                
                my $feature = shift @chr_features;
                if ($feature->{rend} >= $position_lend) {
                    
                    push (@container, $feature);
                }
                else {
                    # no overlap of feature with read range
                    # no op, feature gets tossed.
                }
            }
        
            ## purge current contained features that no longer overlap position.
            @container = sort {$a->{rend}<=>$b->{rend}} @container;
            while (@container && $container[0]->{rend} < $position_lend) {
                shift @container;
            }
            
            
            my $mapped_read_flag = 0;
            if (@container) {
                
                foreach my $feature (@container) {
                    my $acc = $feature->{acc};
                    
                    my $strand_mapping;
                    if ($IGNORE_STRANDEDNESS || $feature->{strand} eq '?') {
                        $strand_mapping = '?'; # for intergenic, consider as sense
                    }
                    elsif ($feature->{strand} eq $trans_feature_strand) {
                        $strand_mapping = "SENSE";
                    }
                    else {
                        $strand_mapping = "ANTI";
                    }

                    
                    if ($NO_REQUIRE_COMPATIBILITY || &Overlap_info::compatible_overlap($feature->{coordset}, $trans_feature_coordset)) {

                        ## determine percent of feature aligning
                        my $overlapping_bases = &Overlap_info::sum_overlaps($feature->{coordset}, $trans_feature_coordset);
                        
                        if ($overlapping_bases > 0) {
                            
                            my $percent_of_trans_length_aligned = sprintf("%.2f", $overlapping_bases / $trans_feature->{length} * 100);
                            
                            print $ofh join("\t", $trans_feature->{acc}, $scaffold, $position_lend, $acc, $strand_mapping, $percent_of_trans_length_aligned) . "\n";
                            
                            $mapped_read_flag = 1;
                            
                        }
                    }

                    
                }
                
            }
            
            
        
            unless ($mapped_read_flag) {
                print $ofh join("\t", $trans_feature->{acc}, ".", ".", ".", ".", ".", ".") . "\n";
            }
            
            #if ($read_counter > 100000) { last; } ###DEBUG
        }
    }
    
    close $ofh;

    return($trans_mapping_file);
}


####
sub refine_mappings_estimate_counts {
    my ($trans_mapping_file) = @_;

    ## sort by read name, molecule, and position
    my $read_name_sorted_file = "$trans_mapping_file.read_name_sorted";
    my $cmd = "sort -k1,1 $trans_mapping_file > $read_name_sorted_file";
    &process_cmd($cmd);


    my @curr_read_structs;
    
    my %feature_to_counts;

    open (my $fh, $read_name_sorted_file) or die "Error, cannot read file $read_name_sorted_file";
    while (<$fh>) {
        chomp;
        my ($read_acc, $mol, $pos, $feature_acc, $sense_or_anti, $percent_mapped) = split(/\t/);
                
        my $read_struct = { 
            read_acc => $read_acc,
            mol => $mol,
            pos => $pos, 
            feature_acc => $feature_acc,
            sense_or_anti => $sense_or_anti,
            percent_mapped => $percent_mapped,
        };
        
        if (@curr_read_structs && $curr_read_structs[0]->{read_acc} ne $read_acc) {
            my @features = &refine_read_mappings(@curr_read_structs);
            foreach my $feature (@features) {
                print join("\t", 
                           $feature->{read_acc},
                           $feature->{mol},
                           $feature->{pos},
                           $feature->{feature_acc},
                           $feature->{sense_or_anti},
                           $feature->{percent_mapped},
                           ) . "\n";
            }
            
            @curr_read_structs = (); # reinit
        }
        
        push (@curr_read_structs, $read_struct);
    }
    close $fh;
    
    if (@curr_read_structs) {
        my @features = &refine_read_mappings(@curr_read_structs);
        
        foreach my $feature (@features) {
            print join("\t", 
                       $feature->{read_acc},
                       $feature->{mol},
                       $feature->{pos},
                       $feature->{feature_acc},
                       $feature->{sense_or_anti},
                       $feature->{percent_mapped},
                       ) . "\n";
            
        }
    }
    
    
    unlink($trans_mapping_file, $read_name_sorted_file); # remove intermediate files
    

}



####
sub refine_read_mappings {
    my (@curr_trans_structs) = @_;
    
    @curr_trans_structs = reverse sort {$a->{percent_mapped}<=>$b->{percent_mapped}} @curr_trans_structs;
    
    if ($IGNORE_ANTISENSE) {
        @curr_trans_structs = grep { $_->{sense_or_anti} !~ /anti/i } @curr_trans_structs;
    }

    unless (@curr_trans_structs) {
        return (); # nothing to do
    }

    if ($BEST_ONLY) {
        @curr_trans_structs = shift @curr_trans_structs;
    }


    return(@curr_trans_structs);


}



####
sub get_coordset_for_isoform {
    my ($isoform) = @_;

    my @coordsets;
    
    my @exons = $isoform->get_exons();
						
    foreach my $exon (@exons) {
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();

        push (@coordsets, [$lend, $rend]);
    }

    @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;


    return(@coordsets);
}


####
sub sum_coordset_segments {
    my ($coordset_aref) = @_;

    my $sum = 0;
    foreach my $coordset (@$coordset_aref) {
        
        my ($lend, $rend) = @$coordset;
        
        $sum += $rend - $lend + 1;
    }

    return($sum);
}

####
sub process_cmd {
    my ($cmd) = @_;
    
    print STDERR "CMD: $cmd\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


####
sub get_min_max_coords {
    my ($read_align_coords_aref) = @_;

    my @coords;
    foreach my $coordset (@$read_align_coords_aref) {
        my ($lend, $rend) = @$coordset;

        push (@coords, $lend, $rend);
    }

    @coords = sort {$a<=>$b} @coords;

    my $min = shift @coords;
    my $max = pop @coords;

    return($min, $max);
}

####
sub compute_aligned_read_length {
    my ($read_coords_aref) = @_;

    my $sum_len = 0;
    foreach my $coordset (@$read_coords_aref) {
        my ($lend, $rend) = @$coordset;
        
        $sum_len += abs($rend - $lend) + 1;
    }

    return($sum_len);
}


####
sub get_feature_lengths {
    my ($chr_to_features_href) = @_;

    my %feature_lengths;

    foreach my $feature_list_aref (values %$chr_to_features_href) {

        foreach my $feature (@$feature_list_aref) {

            my $acc = $feature->{acc};
            my $len = $feature->{length};

            $feature_lengths{$acc} = $len;
        }
    }


    return(%feature_lengths);
}


####
sub merge_short_indels {
    my ($align_coords_aref) = @_;

    if (scalar (@$align_coords_aref) == 1) {
        return($align_coords_aref); # nothing to do
    }

    my @merged_coords = shift @$align_coords_aref;

    foreach my $coordset (shift @$align_coords_aref) {
        my ($lend, $rend) = @$coordset;
        my $prev_rend = $merged_coords[$#merged_coords]->[1];
        if (abs($lend - $prev_rend) <= $MAX_MERGE_INDEL) {
            if ($rend > $prev_rend) {
                $merged_coords[$#merged_coords]->[1] = $rend;
            }
        }
        else {
            # new coord segment
            push (@merged_coords, $coordset);
        }
    }

    return(\@merged_coords);
}


