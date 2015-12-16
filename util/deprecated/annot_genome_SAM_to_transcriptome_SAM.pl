#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../PerlLib/");
use Gene_obj;
use GFF3_utils;
use GTF_utils;
use SAM_reader;
use SAM_entry;
use Overlap_info;
use Fasta_reader;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;
use Storable qw(dclone);
use Nuc_translator;

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
#  --coord_sorted_sam <string>      sam file name
#  
#  --genome <string>                genome file name
#
# # Below are optional:
#  
#  --SS_lib_type       one of [RF,FR,F,R,XS].  if set to XS, then using the XS:A:[+-] attribute for assignment.
#   
#  --out_prefix        output file name (default: 'transcriptome_from_genome' )
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
my $out_prefix = "transcriptome_from_genome";
my $genome_file;

&GetOptions ( 'h' => \$help_flag,
			  'gff3=s' => \$genes_gff3,
			  'gtf=s' => \$genes_gtf,
			  'coord_sorted_sam=s' => \$sam_file,
              'SS_lib_type=s' => \$SS_lib_type,
              'v' => \$VERBOSE,
              'out_prefix=s' => \$out_prefix,
              'genome=s' => \$genome_file,
    );



if ($help_flag) { die $usage; }


unless ($sam_file) {
    die $usage;
}


unless ($genes_gff3 || $genes_gtf) { 
	die $usage;
}

unless ($genome_file && $sam_file) {
    die $usage;
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

    print STDERR "-reading genome sequences\n";
    my $fasta_reader = new Fasta_reader($genome_file);
    my %genome = $fasta_reader->retrieve_all_seqs_hash();
    
    my $cdna_fasta_file = "$out_prefix.cdnas.fasta";
    open (my $cdna_fasta_ofh, ">$cdna_fasta_file") or die $!;
    
    
        
    print STDERR "-organizing transcript features.\n";
    
    ## Populate Genomic Features
	    
	my %chr_to_features;
	
    my %acc_to_seqs;

    my $total_transcripts = 0;

	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
        my $scaffold_seq = $genome{$asmbl_id} or die "Error, cannot find sequence for $asmbl_id";
        

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
                
                my $isoform_id = $isoform->{Model_feat_name};
				my $isoform_length = 0;
				
                my @coordset = &get_coordset_for_isoform($isoform);
                my ($lend, $rend) = ($coordset[0]->[0], $coordset[$#coordset]->[1]);
                
                my $length = &sum_coordset_segments(\@coordset);
                $acc_to_seqs{$isoform_id} = 'N' x $length;
                
                push (@{$chr_to_features{$scaffold}}, { acc => $isoform_id,
                                                        name => $isoform->{com_name},
                                                        lend => $lend,
                                                        rend => $rend,
                                                        coordset => [@coordset],
                                                        strand => $strand,
                                                        length => $length,
                      } );
                
                my $cdna_seq = $isoform->create_cDNA_sequence(\$scaffold_seq);
                
                print $cdna_fasta_ofh ">$isoform_id\n$cdna_seq\n";
                
                

            }
		}
	}
    
    close $cdna_fasta_ofh;
    
	##
	## Map reads to genes
	##
    
	my %gene_id_to_reads_mapped;
	
	my %fragment_tracker;
	my $read_counter = 0;
	my %acc_to_read_count;

	
	## for current processing.
	my $curr_scaff = "";
	my @chr_features;	
	my @container;  ## holds current features.

    
	
    print STDERR "-examining read alignments, mapping to annotated features.\n";
    

    my $trans_mapping_file = "tmp.$out_prefix.pre.unsorted.sam";
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
        my $position_lend = $sam_entry->get_aligned_position();
        
        my ($read_align_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
        
        my $genome_transcribed_orient = "?";
        if ($SS_lib_type) {
            if ($SS_lib_type eq "XS") {
                ## examine XS:A attribute
                my $sam_line = $sam_entry->toString();
                if ($sam_line =~ /XS:A:([\+\-])/) {
                    $genome_transcribed_orient = $1;
                }
            }
            else {
                $genome_transcribed_orient = $sam_entry->get_query_transcribed_strand($SS_lib_type);
            }
        }
        
        
        if ($molecule ne $curr_scaff) {
            ## reset
            if (exists $chr_to_features{$molecule}) {
                
                @chr_features = sort {$a->{lend}<=>$b->{lend}} @{$chr_to_features{$molecule}};
            }
            else {
                print STDERR "Warning, no chr features for molecule: $molecule\n";
            }
            @container = (); # clear
            $curr_scaff = $molecule;
        }
        
        @container = sort {$a->{rend}<=>$b->{rend}} @container;
        ## purge current contained features that no longer overlap position.
        while (@container && $container[0]->{rend} < $position_lend) {
            shift @container;
        }
        
        ## collect features that overlap 
        while (@chr_features && 
               $chr_features[0]->{lend} <= $position_lend) {
            
            my $feature = shift @chr_features;
            if ($feature->{rend} >= $position_lend) {
                
                push (@container, $feature);
            }
            else {
                # no op, feature gets tossed.
            }
        }
        
        if (@container) {
            
            foreach my $feature (@container) {
                my $acc = $feature->{acc};
                if ( #($FUZZY_OVERLAP && &Overlap_info::compatible_overlap($feature->{coordset}, $read_align_coords_aref) )
                     # ||
                      (&Overlap_info::compatible_overlap_A_contains_B($feature->{coordset}, $read_align_coords_aref))
                      ) {
                    
                    print STDERR "\r[$read_counter] reads examined.     Adding read count to $acc  " if $VERBOSE;
                    
                    &write_transcript_sam_record($ofh, $sam_entry, $genome_transcribed_orient, $feature, $read_align_coords_aref);
                    
                }
                
            }
            
            
            
            
        }
        else {
            #print STDERR "\nNo annot mapped to $_\n";
            
        }
        
        
        #if ($read_counter > 100000) { last; } ###DEBUG
    }
    
    
    close $ofh;
    
    ## sort by read name, scaff name and coordinate
    my $coord_sorted_pre_sam = "$out_prefix.pre.nameSorted.sam";
    my $cmd = "sort -k1,1 -k3,3 -k4,4n $trans_mapping_file > $coord_sorted_pre_sam";
    &process_cmd($cmd);
    
    my $refined_pairing_file = "$out_prefix.nameSorted.sam";
    &reestablish_read_pairing($coord_sorted_pre_sam, $refined_pairing_file);
    

    ## delete the intermediates
    
    #unlink($trans_mapping_file);
    #unlink($coord_sorted_pre_sam);


    # create bam file for use with RSEM
    
    $cmd = "samtools faidx $cdna_fasta_file";
    &process_cmd($cmd);

    $cmd = "samtools view -bt $cdna_fasta_file.fai $refined_pairing_file > $out_prefix.nameSorted.bam";
    &process_cmd($cmd);
    
        
    print "\n\nFinished.  See file: $refined_pairing_file\n\n";
    

    exit(0);
}


####
sub generate_expressed_transcript_report {
    my ($fpkm_file, $total_transcripts) = @_;

    my $Rscript = "$fpkm_file.R";
    open (my $ofh, ">$Rscript");
    print $ofh "source(\"$FindBin::RealBin/R/expression_analysis_lib.R\")\n";
    print $ofh "png(\"$fpkm_file.genes_vs_minFPKM.png\", width=1000, height=500)\n";
    print $ofh "plot_expressed_gene_counts(\"$fpkm_file\", total=$total_transcripts, title=\"expressed transcripts vs. min FPKM\", fpkm_range=seq(0,5,0.01), outfile=\"$fpkm_file.genes_vs_minFPKM.dat\")\n";
    print $ofh "dev.off()\n";
    close $ofh;
    
    system("R --vanilla -q < $Rscript");
    
    return;
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
sub write_transcript_sam_record {
    my ($ofh, $sam_entry, $genome_transcribed_orient, $feature, $read_align_coords_aref) = @_;

    $sam_entry = dclone($sam_entry); # make a deep copy, don't change the original record.
    
    #print Dumper($feature);
    #print Dumper($read_align_coords_aref);

    my $cdna_coords_aref = &make_cdna_coords($read_align_coords_aref, $feature);
    
    ## KISS here too, single segment for now, ignoring gaps
    my $left_coord = $cdna_coords_aref->[0];
    
    if ($feature->{strand} eq '-') {
        # must swap the read orientation
        my $query_strand = $sam_entry->get_query_strand();
        my $opposite_strand = ($query_strand eq '+') ? '-' : '+';
        $sam_entry->set_query_strand($opposite_strand);
    }

    my @fields = $sam_entry->get_fields();
    
    $fields[2] = $feature->{acc};
    $fields[3] = $left_coord;
    
    my $read_length = length($sam_entry->get_sequence());
    $fields[5] = "${read_length}M"; 
    
    if ($feature->{strand} eq '-') {
        
        # revcomp the sequence and the quals

        $fields[9] = &reverse_complement($fields[9]);
        $fields[10] = join("", reverse( split(//, $fields[10]) ) );

    }
    

    my $right_pos = $left_coord + $read_length - 1;
    #print "$left_coord-$right_pos of " . $feature->{length} . "\n";
    
    if ($left_coord + $read_length - 1 > $feature->{length}) {
        
        die "Error, read: $left_coord-$right_pos exceeds transcript length: " . $feature->{length} . "  ";
    }

    
    # remove XS field if exists
    @fields = grep { $_ !~ /^XS:A:/ } @fields;
    
    
    if ($genome_transcribed_orient ne '?') {
       
        my $transcribed_orient = ($genome_transcribed_orient);
        if ($feature->{strand} eq '-') {
            # flip it
            $transcribed_orient = ($transcribed_orient eq '+') ? '-' : '+';
        }
        
        push (@fields, "XS:A:$transcribed_orient");
    }

    print $ofh join("\t", @fields) . "\n";
    
    
    return;
    
}


####
sub make_cdna_coords {
    my ($read_coords_aref, $feature) = @_;

    ## TODO: full proper coordinate transformation. Right now, KISS and use boundaries
    
    my $read_lend = $read_coords_aref->[0]->[0];
    my $read_rend = $read_coords_aref->[$#$read_coords_aref]->[1];

    my $trans_lend_coord = &convert_genome_to_cdna_coord($read_lend, $feature->{coordset});
    my $trans_rend_coord = &convert_genome_to_cdna_coord($read_rend, $feature->{coordset});
    
    if ($feature->{strand} eq '-') {
        # revcomp the coords
        $trans_lend_coord = $feature->{length} - $trans_rend_coord + 1;
        $trans_rend_coord = $feature->{length} - $trans_lend_coord + 1;
    }

    return([$trans_lend_coord, $trans_rend_coord]);
}


####
sub convert_genome_to_cdna_coord {
    my ($coord_to_convert, $genome_coords_aref) = @_;

    my $ret_coord = 0;

    foreach my $coordset (@$genome_coords_aref) {
        
        my ($lend, $rend) = @$coordset;
        
        if ($coord_to_convert >= $lend && $coord_to_convert <= $rend) {

            $ret_coord += $coord_to_convert - $lend + 1;
            
            return($ret_coord);
        }
        else {
            $ret_coord += $rend - $lend + 1; # just the segment length 
        }
    }

    die "Error, couldn't convert coordinate: $coord_to_convert within " . Dumper($genome_coords_aref);
}

####
sub reestablish_read_pairing {
    my ($input_sam, $output_sam) = @_;
    
    open (my $ofh, ">$output_sam") or die "Error, cannot write to $output_sam";
    
    my $sam_reader = new SAM_reader($input_sam);
    
    my $sam_entry_A = $sam_reader->get_next();
    my $sam_entry_B = $sam_reader->get_next();

    while ($sam_entry_A) {

        if (&same_read_name($sam_entry_A, $sam_entry_B) && &opposite_fragment_ends($sam_entry_A, $sam_entry_B)) {
            &make_proper_pairs($sam_entry_A, $sam_entry_B);
            print $ofh $sam_entry_A->toString() . "\n";
            print $ofh $sam_entry_B->toString() . "\n";
            $sam_entry_A = $sam_reader->get_next();
            $sam_entry_B = $sam_reader->get_next();
            
        }
        else {
            &make_unpaired($sam_entry_A);
            print $ofh $sam_entry_A->toString() . "\n";
            $sam_entry_A = $sam_entry_B;
            $sam_entry_B = $sam_reader->get_next();
        }
    }

    close $ofh;

    return;
}
        
####
sub same_read_name {
    my ($sam_entry_A, $sam_entry_B) = @_;

    if (! ($sam_entry_A && $sam_entry_B) ) {
        return(0);
    }

    if ($sam_entry_A->get_read_name() eq $sam_entry_B->get_read_name()) {
        return(1);
    }
    else {
        return(0);
    }
}

####
sub opposite_fragment_ends {
    my ($sam_entry_A, $sam_entry_B) = @_;

    if (! ($sam_entry_A && $sam_entry_B) ) {
        return(0);
    }


    if (  ($sam_entry_A->is_first_in_pair() && $sam_entry_B->is_second_in_pair())
          ||
          ($sam_entry_A->is_second_in_pair() && $sam_entry_B->is_first_in_pair())   ) {

        return(1);
    }
    else {
        return(0);
    }
}

####
sub make_proper_pairs {
    my ($sam_entry_A, $sam_entry_B) = @_;

    $sam_entry_A->set_proper_pair(1);
    $sam_entry_B->set_proper_pair(1);
    
    $sam_entry_A->set_mate_scaffold_name($sam_entry_B->get_scaffold_name());
    $sam_entry_A->set_mate_scaffold_position($sam_entry_B->get_scaffold_position());

    $sam_entry_B->set_mate_scaffold_name($sam_entry_A->get_scaffold_name());
    $sam_entry_B->set_mate_scaffold_position($sam_entry_A->get_scaffold_position());
    
    $sam_entry_A->set_mate_unmapped(0);
    $sam_entry_B->set_mate_unmapped(0);

    
    return;
}

####
sub make_unpaired {
    my ($sam_entry) = @_;

    $sam_entry->set_proper_pair(0);
    if ($sam_entry->is_paired()) {
        $sam_entry->set_mate_unmapped(1);
    }
    
    $sam_entry->set_mate_scaffold_name("*");
    $sam_entry->get_mate_scaffold_position("*");
    
    return;
}
