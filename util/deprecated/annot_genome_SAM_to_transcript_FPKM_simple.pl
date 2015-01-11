#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib/");
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
#  ### Annot settings:
#
#  --gff3 <string>     gff3 file name
# 
#    OR
#
#  --gtf <string>      gtf file name 
#
#    AND
#
# ### SAM settings:
#
#  --coord_sorted_sam <string>  genome coordinate-sorted sam file name
#
# ### rna-seq info
#
#  --frag_length      length of an RNA-Seq fragment (** not the length of a read, but rather the mean length of a fragment from which the read was derived **)
#
# # Below are optional:
#  
#  --SS_lib_type       one of [RF,FR,F,R,XS].  if set to XS, then using the XS:A:[+-] attribute for assignment.
#   
#  --FUZZY             require compatible overlap of read/transcript but do not require containment of entire read.
#  --REQUIRE_PAIRINGS  transcript must have both paired fragment ends mapped in order to be counted.
#
#
#  --outfile           output file name (default:  sam_filename.fpkm)
#
#  --include_intergenic   create features out of the intergenic regions.
#  --extend_UTRs <int>    default: 0
#
###################################################################################################


_EOUSAGE_

	;




my $genes_gff3;
my $genes_gtf;
my $sam_file;
my $help_flag;
my $SS_lib_type;
my $VERBOSE = 0;
my $FUZZY_OVERLAP = 0;
my $REQUIRE_PAIRINGS = 0;
my $outfile = "";
my $frag_length;;
my $INCLUDE_INTERGENIC = 0;

my $MIN_INTERGENIC_LENGTH = 100;
my $MIN_EFF_LEN = 10;
my $MAX_MERGE_INDEL = 5;
my $EXTEND_UTR_LENGTH = 0;


&GetOptions ( 'h' => \$help_flag,
			  'gff3=s' => \$genes_gff3,
			  'gtf=s' => \$genes_gtf,
			  'coord_sorted_sam=s' => \$sam_file,
              'SS_lib_type=s' => \$SS_lib_type,
              'v' => \$VERBOSE,
              'FUZZY' => \$FUZZY_OVERLAP,
              'REQUIRE_PAIRINGS' => \$REQUIRE_PAIRINGS,
              'outfile=s' => \$outfile,
              'frag_length=i' => \$frag_length,
              
              'include_intergenic' => \$INCLUDE_INTERGENIC,
              
              'extend_UTRs=i' => \$EXTEND_UTR_LENGTH,

              );

if (@ARGV) {
    die "Error, cannot parse param: @ARGV";
}


if ($help_flag) { die $usage; }

unless ($sam_file) {
    die $usage;
}

unless ($genes_gff3 || $genes_gtf) { 
	die $usage;
}


unless ($frag_length) {
    die $usage;
}

unless ($outfile) {
    $outfile = "$sam_file.fpkm";
}


main: {
	
    my $gene_obj_indexer_href = {};

	my $contig_to_gene_list_href;

    print STDERR "-parsing gene annotations file\n";
	if ($genes_gff3) {
	
		## associate gene identifiers with contig id's.
		$contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($genes_gff3, $gene_obj_indexer_href);
	}
	else {
		# GTF mode
		$contig_to_gene_list_href = &GTF_utils::index_GTF_gene_objs_from_GTF($genes_gtf, $gene_obj_indexer_href);
	}
    
    
  
    
    ## Populate Genomic Features

    print STDERR "-organizing annotated transcript features.\n";

	my %chr_to_features = &get_chr_to_features($contig_to_gene_list_href, $gene_obj_indexer_href);

    if ($INCLUDE_INTERGENIC) {
        print STDERR "-adding intergenic features\n";
        &add_intergenic_features(\%chr_to_features);
    }

	my %feature_lengths = &get_feature_lengths(\%chr_to_features);
    
    
    ## assign read mapping to features, including sense and antisense mappings
    print STDERR "\n\n-mapping reads to features\n";
    my $mapped_read_file = &map_reads_to_features(\%chr_to_features, $sam_file);

    print STDERR "\n\n-refining mappings and estimating counts\n\n";
    my ($total_mapped_fragments, $feature_to_counts_href, $refined_mapping_file) = &refine_mappings_estimate_counts($mapped_read_file); # count fragments, not individual reads
    
    print STDERR "\n\n-generating summary output\n";
    &report_summary_file($feature_to_counts_href, $total_mapped_fragments, \%feature_lengths, $frag_length, $outfile);
    

    exit(0);
    

}


####
sub get_chr_to_features {
    my ($contig_to_gene_list_href, $gene_obj_indexer_href) = @_;


    my %chr_to_features;
	
    my $total_transcripts = 0;

	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};

        my @exon_coords_plus;
        my @exon_coords_minus;

		foreach my $gene_id (@gene_ids) {
			my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
			
			unless (ref $gene_obj_ref) {
				die "Error, no gene_obj for gene_id: $gene_id";
			}
			
			my $strand = $gene_obj_ref->get_orientation();

			
			my $scaffold = $asmbl_id;
            
			my $max_isoform_cdna_length = 0;
			
            my $exon_coords_aref = ($strand eq '+') ? \@exon_coords_plus : \@exon_coords_minus;
            
            
			foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                
                $total_transcripts++;
                
                my $isoform_id = join("::", $gene_id, $isoform->{Model_feat_name}); ## embed the gene identifier in with the transcript id, so we can tease it apart later.
				my $isoform_length = 0;
				
                my @coordset = &get_coordset_for_isoform($isoform);
                push (@$exon_coords_aref, @coordset);
                
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
                        acc => "intergenic_region.$inter_counter",
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
sub map_reads_to_features {
    my ($chr_to_features_href, $sam_file) = @_;
    
	##
	## Map reads to genes
	##

	my %gene_id_to_reads_mapped;
	
	my $read_counter = 0;
	
	## for current processing.
	my $curr_scaff = "";
	my @chr_features;	
	my @container;  ## holds current features.

    
	
    print STDERR "-examining read alignments, mapping to annotated features.\n";
    

    my $trans_mapping_file = "read_trans_mappings.$$.txt";
    open (my $ofh, ">$trans_mapping_file") or die $!;
    
    
    my $sam_reader = new SAM_reader($sam_file);
    while ($sam_reader->has_next()) {
        
        my $sam_entry = $sam_reader->get_next();
                
        if ($sam_entry->is_query_unmapped()) { 
            next;
        }
        
        $read_counter++;
        
        if ($read_counter % 1000 == 0) {
            print STDERR "\r[$read_counter] reads examined.     ";
        }
        
        my $read_acc =  $sam_entry->get_core_read_name();
        my $full_read_acc = $sam_entry->reconstruct_full_read_name();
                
        my $molecule = $sam_entry->get_scaffold_name();
        #my $position_lend = $sam_entry->get_aligned_position();
        
        my ($read_align_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
        
        $read_align_coords_aref = &merge_short_indels($read_align_coords_aref);
        
        my ($position_lend, $position_rend) = &get_min_max_coords($read_align_coords_aref);
        my $aligned_read_length = &compute_aligned_read_length($read_align_coords_aref);
        
        my $orient = "?";
        if ($SS_lib_type) {
            if ($SS_lib_type eq "XS") {
                ## examine XS:A attribute
                my $sam_line = $sam_entry->toString();
                if ($sam_line =~ /XS:A:([\+\-])/) {
                    $orient = $1;
                }
            }
            else {
                $orient = $sam_entry->get_query_transcribed_strand($SS_lib_type);
            }
        }
        
        
              			
        
        if ($molecule ne $curr_scaff) {
            ## reset
            if (exists $chr_to_features_href->{$molecule}) {
                
                @chr_features = sort {$a->{lend}<=>$b->{lend}} @{$chr_to_features_href->{$molecule}};
            }
            else {
                print STDERR "Warning, no chr features for molecule: $molecule\n";
            }
            @container = (); # clear
            $curr_scaff = $molecule;
        }
        

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
                    
                my $strand_mapping = ($orient eq "?" || $feature->{strand} eq $orient) ? "S" : "A";
                
            
                if (&Overlap_info::compatible_overlap_A_contains_B($feature->{coordset}, $read_align_coords_aref)) {
                    
                    print $ofh join("\t", $full_read_acc, $molecule, $position_lend, $acc, "C", $strand_mapping) . "\n";
                    $mapped_read_flag = 1;
                    
                }
                
                elsif ($FUZZY_OVERLAP 
                       && 
                       &Overlap_info::compatible_overlap($feature->{coordset}, $read_align_coords_aref) 
                       &&
                       &Overlap_info::sum_overlaps($feature->{coordset}, $read_align_coords_aref) >= $aligned_read_length/2) {
                                        
                    print $ofh join("\t", $full_read_acc, $molecule, $position_lend, $acc, "P", $strand_mapping) . "\n";
                    $mapped_read_flag = 1;
                }
                elsif ($strand_mapping eq "A" # antisense
                       &&
                       &Overlap_info::sum_overlaps($feature->{coordset}, $read_align_coords_aref) >= $aligned_read_length/2) {
                    
                    ## antisense mapping less strict than sense-mapping, compatibility not required.
                    ## more than half the read aligns as antisense to current feature, but otherwise incompatible
                    print $ofh join("\t", $full_read_acc, $molecule, $position_lend, $acc, "P", $strand_mapping) . "\n";
                }
                
            }
                        
        }
        
        unless ($mapped_read_flag) {
            print $ofh join("\t", $full_read_acc, ".", ".", ".", ".", ".") . "\n";
        }
        
        #if ($read_counter > 100000) { last; } ###DEBUG
    }
    
    
    close $ofh;

    return($trans_mapping_file);
}


####
sub refine_mappings_estimate_counts {
    my ($trans_mapping_file) = @_;

    ## sort by read name, molecule, and position
    my $read_name_sorted_file = "$trans_mapping_file.read_name_sorted";
    my $cmd = "sort -k1,1 -k2,2 -k3,3n $trans_mapping_file > $read_name_sorted_file";
    &process_cmd($cmd);


    my $count_reads = 0;

    my $refined_file = "$read_name_sorted_file.refined";
    open (my $ofh, ">$refined_file") or die "Error, cannot write to file $refined_file";
    

    my @curr_read_structs;
    
    my %feature_to_counts;

    open (my $fh, $read_name_sorted_file) or die "Error, cannot read file $read_name_sorted_file";
    while (<$fh>) {
        chomp;
        my ($read_acc, $mol, $pos, $feature_acc, $complete_or_partial, $sense_or_anti) = split(/\t/);
        
        my $core_read_acc = $read_acc;
        $core_read_acc =~ s|/\d+$||;

        if ($sense_or_anti eq "A") {
            $feature_acc .= "__ANTI";  ## very hacky... need to change  TODO
        }
        
        my $read_struct = { 
            core_read_acc => $core_read_acc,
            read_acc => $read_acc,
            mol => $mol,
            pos => $pos, 
            feature_acc => $feature_acc,
            sense_or_anti => $sense_or_anti,
            complete_or_partial => $complete_or_partial,
        };
        
        if (@curr_read_structs && $curr_read_structs[0]->{core_read_acc} ne $core_read_acc) {
            my @features = &refine_read_mappings($ofh, \@curr_read_structs);
            $count_reads++;
            foreach my $feature (@features) {
                #print STDERR "incrementing $feature\n";
                $feature_to_counts{$feature} += 1/(scalar @features);
            
                #print Dumper(\%feature_to_counts);
            }
            
            @curr_read_structs = (); # reinit
            
            if ($count_reads % 1000 == 0) {
                print STDERR "\r[$count_reads] frags processed   ";
            }
        }
        
        push (@curr_read_structs, $read_struct);
    }
    close $fh;

    if (@curr_read_structs) {
        my @features = &refine_read_mappings($ofh, \@curr_read_structs);
        $count_reads++;
        foreach my $feature (@features) {
            $feature_to_counts{$feature} += 1/(scalar @features);
        }
    }
    
    close $ofh;
    
    return($count_reads, \%feature_to_counts, $refined_file);
}




####
sub refine_read_mappings {
    my ($ofh, $curr_read_structs_aref) = @_;


    #print STDERR "\rrefining mappings for: " . $curr_read_structs_aref->[0]->{core_read_acc} . "   ";
    
    ## group reads by feature

    my %feature_to_struct;

    foreach my $struct (@$curr_read_structs_aref) {
        my $feature_acc = $struct->{feature_acc};
        if ($feature_acc ne ".") {
            push (@{$feature_to_struct{$feature_acc}}, $struct);
        }
    }

    my @paired_structs;
    my @unpaired_structs;
    foreach my $struct_set_aref (values %feature_to_struct) {
        if (scalar(@$struct_set_aref) == 2) {
            push (@paired_structs, @$struct_set_aref);
        }
        else {
            push (@unpaired_structs, @$struct_set_aref);
        }
    }
    
    ## report all paired feature data
    my %counted_features;


    my %mol_pos_seen;
    foreach my $paired_struct (@paired_structs) {
        print $ofh join("\t", $paired_struct->{read_acc}, $paired_struct->{mol}, $paired_struct->{pos}, 
                        $paired_struct->{feature_acc}, $paired_struct->{sense_or_anti}, $paired_struct->{complete_or_partial}) . "\n";
    
        my $mol_pos_key = join(";", $paired_struct->{mol}, $paired_struct->{pos});
        $mol_pos_seen{$mol_pos_key}=1;
        
        $counted_features{ $paired_struct->{feature_acc} }++;
    
    }


    if ($REQUIRE_PAIRINGS) {
        ## done here. Ignore any of the reads that do not map as pairs to a given feature
        return (keys %counted_features);
    }


    ## prioritize unpaired by complete status.  
    my @complete_unpaired = grep {$_->{complete_or_partial} eq "C"} @unpaired_structs;

    foreach my $complete_unpaired (@complete_unpaired) {

        my $mol_pos_key = join(";", $complete_unpaired->{mol}, $complete_unpaired->{pos});
        $mol_pos_seen{$mol_pos_key}=1;


        print $ofh join("\t", $complete_unpaired->{read_acc}, $complete_unpaired->{mol}, $complete_unpaired->{pos}, 
                        $complete_unpaired->{feature_acc}, $complete_unpaired->{sense_or_anti}, $complete_unpaired->{complete_or_partial}) . "\n";
    

        $counted_features{ $complete_unpaired->{feature_acc} }++;

    }

    ## report remaining mol/pos not paired
    my @partial_unpaired = grep {$_->{complete_or_partial} eq "P"} @unpaired_structs;
    foreach my $unpaired (@partial_unpaired) {
        my $mol_pos_key = join(";", $unpaired->{mol}, $unpaired->{pos});
        if ($mol_pos_seen{$mol_pos_key}) {
            ## already mapped to a feature as complete and/or as paired
            next;
        }
        
        print $ofh join("\t", $unpaired->{read_acc}, $unpaired->{mol}, $unpaired->{pos}, 
                        $unpaired->{feature_acc}, $unpaired->{sense_or_anti}, $unpaired->{complete_or_partial}) . "\n";
        
        

        $counted_features{ $unpaired->{feature_acc} }++;

    }
    
    return (keys %counted_features);
}



####
sub get_coordset_for_isoform {
    my ($isoform) = @_;

    my @coordsets;
    
    my @exons = sort {$a->{end5}<=>$b->{end5}} $isoform->get_exons();
    
    for (my $i = 0; $i <= $#exons; $i++) { 
        
        my $exon = $exons[$i];
        my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
        
        if ($EXTEND_UTR_LENGTH && $i == 0) {
            $lend -= $EXTEND_UTR_LENGTH;
        }
        if ($EXTEND_UTR_LENGTH && $i == $#exons) {
            $rend += $EXTEND_UTR_LENGTH;
        }
        
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
sub report_summary_file {
    my ($feature_to_counts_href, $total_mapped_fragments, $feature_lengths_href, $fragment_length, $outfile) = @_;


    open (my $ofh, ">$outfile") or die "Error, cannot open file $outfile for writing";
    
    my $report_text = "";
    $report_text .= "## Total reads mapped: $total_mapped_fragments\n";
    $report_text .= join("\t", "#acc", "count", "length", "eff_length", "FPKM") . "\n";
    
    foreach my $feature (sort keys %$feature_to_counts_href) {
        my $count = $feature_to_counts_href->{$feature};
        
        my $feat_name = $feature;
        if ($feat_name =~ /__ANTI$/) {
            $feat_name =~ s/__ANTI$//;
        }
        
        my $length = $feature_lengths_href->{$feat_name};
        
        my $eff_length = $length - $fragment_length + 1;
        if ($eff_length < $MIN_EFF_LEN) {
            $eff_length = $MIN_EFF_LEN;
        }


        my $fpkm = sprintf("%.2f", $count / ($length/1e3) / ($total_mapped_fragments/1e6));


        $report_text .= join("\t", $feature, $count, $length, $eff_length, $fpkm) . "\n";
        
    }
    
    print $report_text;
    print $ofh $report_text;

    close $ofh;
    

    print "\n\nDone.  See file: $outfile\n";

    return;
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

