#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../../PerlLib");
use Gene_obj;
use GFF3_utils;
use BED_utils;
use Carp;
use Nuc_translator;
use Fasta_reader;

my $usage = "\n\nusage: $0 reference.gff3,bed target.gff3,bed [genome.fa restricts to introns w/ consensus splices]\n\n";

my $reference_gff3_bed = $ARGV[0] or die $usage;
my $target_gff3_bed = $ARGV[1] or die $usage;
my $genome_fa = $ARGV[2];

my $DEBUG = 0;

unless ($reference_gff3_bed =~ /(\.bed|\.gff3)$/) {
	die $usage;
}

unless ($target_gff3_bed =~ /(\.bed|\.gff3)$/) {
    die $usage;                                                                                                                               
}                                                   


my %genome;
if ($genome_fa) {

	my $fasta_reader = new Fasta_reader($genome_fa);
	%genome = $fasta_reader->retrieve_all_seqs_hash();
}


my %reference_introns;
my %reference_combos;
my %reference_intron_to_combos;

{ # parse reference:

	my $gene_obj_indexer_href = {};
	
	my $contig_to_gene_list_href = ($reference_gff3_bed =~ /\.gff3$/) 
		? &GFF3_utils::index_GFF3_gene_objs($reference_gff3_bed, $gene_obj_indexer_href)
		: &BED_utils::index_BED_as_gene_objs($reference_gff3_bed, $gene_obj_indexer_href);
	
	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
				
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
		
		foreach my $gene_id (@gene_ids) {
			my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
					
			foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
 
				my $orientation = $isoform->get_orientation();
					
				my $gene_id = $isoform->{TU_feat_name};
				my $isoform_id = $isoform->{Model_feat_name};
				
				my @introns = $isoform->get_intron_coordinates();
				@introns = sort {$a->[0]<=>$b->[0]} @introns;
				
				if (@introns) {
					my $intron_text = "$asmbl_id:"; 
					
					my @intron_labels;
					foreach my $intron (@introns) {
						my ($end5, $end3) = sort {$a<=>$b} @$intron;
						
						$intron_text .= "_$end5-${end3}_";
						
						my $intron_label = "$asmbl_id:$end5-${end3}";

						$reference_introns{$intron_label} .= "$gene_id,$isoform_id;";
						push (@intron_labels, $intron_label);
					}
					
					$reference_combos{$intron_text} .= "$gene_id,$isoform_id;";
					foreach my $intron_label (@intron_labels) {
						$reference_intron_to_combos{$intron_label}->{$intron_text} = 1;
					}
					
				}
			}
		}
	}
}


if ($DEBUG) {
	open (my $ofh, ">ref_info.txt") or die $!;
	foreach my $key (keys %reference_combos) {
		my $val = $reference_combos{$key};
		print $ofh "$key\t$val\n";
	}

	foreach my $key (keys %reference_introns) {
		my $val  = $reference_introns{$key};
		print $ofh "$key\t$val\n";
	}
}




## examine target:

my $gene_obj_indexer_href = {};

my $contig_to_gene_list_href = ($target_gff3_bed =~ /\.gff3$/) 
	? &GFF3_utils::index_GFF3_gene_objs($target_gff3_bed, $gene_obj_indexer_href)
	: &BED_utils::index_BED_as_gene_objs($target_gff3_bed, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
	
	my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};


	my $genome_seq;
	if (%genome) {
		$genome_seq = $genome{$asmbl_id};
		unless ($genome_seq) {
			die "Error, cannot find genome sequence for [$asmbl_id]";
		}
	}
	
	my $genome_seq_sref = \$genome_seq;
	
	foreach my $gene_id (@gene_ids) {
		my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
		my $orient = $gene_obj_ref->get_orientation();
		
		foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
			
			my $orientation = $isoform->get_orientation();
			
			my $gene_id = $isoform->{TU_feat_name};
			my $isoform_id = $isoform->{Model_feat_name};
			
			my @introns = $isoform->get_intron_coordinates();
			
			@introns = sort {$a->[0]<=>$b->[0]} @introns;

			if (@introns) {
				
				my $INTRONS_VALID_FLAG = 1;

				my @intron_labels;
				
				my $intron_text = "$asmbl_id:"; 
				foreach my $intron (@introns) {
					my ($end5, $end3) = sort {$a<=>$b} @$intron;
					
					#my $intron_seq = substr($$genome_seq_sref, $end5-1, $end3-$end5+1);
					#print "intron ($end5-$end3,$orient): $intron_seq\n";

					$intron_text .= "_$end5-${end3}_";
					
					my $intron_label = "$asmbl_id:$end5-${end3}";
					
					if ($genome_seq) {
						unless ($reference_introns{$intron_label} || 
								&consensus_intron_boundaries($genome_seq_sref, $orient, $end5, $end3)) {
							$INTRONS_VALID_FLAG = 0;
							next;
						}
					}
					
					push (@intron_labels, $intron_label);
					
					if (my $gene_info =  $reference_introns{$intron_label}) {
						print "REF\t$asmbl_id:$end5-$end3\t$gene_info\n";
					}
					else {
						print "NOVEL\t$asmbl_id:$end5-$end3\t$gene_id,$isoform_id\n";
					}
				}
				
				if ($INTRONS_VALID_FLAG && scalar @introns > 1) {
					
					if (my $gene_info = $reference_combos{$intron_text}) {
						print "REFCOMBO\t$intron_text\t$gene_info\n";
					}
					else {
						## see if it's a substring of a reference path containing that intron.
						my $found_as_substring = 0;
					  intron_subsearch: 
						foreach my $intron_label (@intron_labels) {
							foreach my $intron_combo (keys %{$reference_intron_to_combos{$intron_label}}) {
								if ($intron_combo =~ /$intron_text/) {
									$found_as_substring = 1;
									last intron_subsearch;
								}
							}
						}
						if ($found_as_substring) {
							print "SUBrefCombo\t$intron_text\t$gene_id,$isoform_id\n";
						}
						else {
							print "NOVELCOMBO\t$intron_text\t$gene_id,$isoform_id\n";
						}
						
					}
					
				}
			}
		}
	}
}



exit(0);

####
sub consensus_intron_boundaries {
	my ($genome_seq_sref, $orient, $lend, $rend) = @_;

	my $seq_len = length($$genome_seq_sref);
	unless ($lend < $seq_len && $rend < $seq_len) {
		print STDERR "-error, coordinates ($lend, $rend) are not within range of sequence: 1-$seq_len\n";
		return (0);
	}

	my $left_splice_boundary = &get_left_splice($genome_seq_sref, $lend);
	my $right_splice_boundary = &get_right_splice($genome_seq_sref, $rend);
	
	#print "SPLICE: $left_splice_boundary, $right_splice_boundary\n";
	
	if ($orient eq '-') {
		($left_splice_boundary, $right_splice_boundary) = ($right_splice_boundary, $left_splice_boundary);
		
		$left_splice_boundary = &reverse_complement($left_splice_boundary);
		$right_splice_boundary = &reverse_complement($right_splice_boundary);

	}

	if ($left_splice_boundary =~ /^(GT|GC)$/ && $right_splice_boundary eq "AG") {
		return(1);
	}
	else {
		return(0);
	}
}

####
sub get_left_splice {
	my ($genome_seq_sref, $coord) = @_;

	my $subseq = substr($$genome_seq_sref, $coord-1, 2);

	return(uc $subseq);
}

####
sub get_right_splice {
	my ($genome_seq_sref, $coord) = @_;

	my $subseq = substr($$genome_seq_sref, $coord-2, 2);

	return(uc $subseq);
}
