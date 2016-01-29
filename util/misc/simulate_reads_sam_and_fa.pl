#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib "$ENV{TRINITY_HOME}/PerlLib";
use Simulate::Uniform_Read_Generator;
use Overlap_info;
use GFF3_utils;
use GTF_utils;
use Gene_obj;
use Fasta_reader;
use SAM_entry;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use CIGAR;
use Nuc_translator;
use List::Util qw (shuffle);
use Storable qw (dclone);

$ENV{LC_ALL} = 'C';


my $usage = <<__EOUSAGE__;

#####################################################################################
#
# Required:
#
#  --gff3 <string>                   reference annotations in gff3 format
#    or
#  --gtf <string>                    reference annotations in gtf format
#
#  --genome <string>                 reference genome in fasta file format
#
#  --frag_length <int>               mean fragment length 
#
#  --read_length <int>               read length from end of fragment
#
#
# Optional:
#
#  --SS_lib_type <string>            if paired: RF, FR;  if single: R or F
#
#  --paired                          default (off, meaning single)
#
#  --out_prefix <string>             default: 'simul'
#
####################################################################################



__EOUSAGE__

    ;

#  --frag_length_stdev <int>         standard deviation for normal distribution of fragment lengths


my $gff3_file;
my $gtf_file;
my $read_length;
my $mean_frag_length;
my $frag_length_stdev = 10; # not used yet

my $reads_per_transcript_file;
my $genome_file;
my $SS_lib_type;
my $paired_flag = 0;
my $out_prefix = "simul";


&GetOptions( "gff3=s" => \$gff3_file,
             "gtf=s" => \$gtf_file,
             "read_length=i" => \$read_length,
             "frag_length=i" => \$mean_frag_length,
             "frag_length_stdev=i" => \$frag_length_stdev,
             
             "genome=s" => \$genome_file,
             "SS_lib_type=s" => \$SS_lib_type,
             "paired" => \$paired_flag,
             "out_prefix=s" => \$out_prefix,
             
			 );

unless ( ($gff3_file || $gtf_file) && $genome_file
         && $mean_frag_length && $read_length && $frag_length_stdev) {
    die $usage;
}


if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(RF|FR|F|R)$/) {
        die "Error, do not recognize SS_lib_type: [$SS_lib_type] ";
    }
    if ($paired_flag) {
        unless ($SS_lib_type eq "RF" || $SS_lib_type eq "FR") {
            die "Error, SS_lib_type: $SS_lib_type is not acceptable for paired reads.";
        }
    }

}


if ($read_length > $mean_frag_length) {
    die "Error, read length: $read_length exceeds the fragment lenth: $mean_frag_length ";
}

main: {

    srand();
    
    my $genome_sam_outfile = "$out_prefix.genome.sam"; 
    open (my $genome_sam_FH, ">$genome_sam_outfile") or die "Error, cannot write to $genome_sam_outfile";    
    
    my $trans_sam_outfile = "$out_prefix.transcriptome.sam";
    open (my $trans_sam_FH, ">$trans_sam_outfile") or die "Error, cannot write to $trans_sam_outfile";
    
    
    my $cdna_file = "$out_prefix.transcriptome.cdnas";
    open (my $cdna_ofh, ">$cdna_file") or die $!;
   
    
    my $fasta_reader = new Fasta_reader($genome_file);
    my %transcript_seqs = $fasta_reader->retrieve_all_seqs_hash();
    
    my $reads_file = "$out_prefix.reads.fa";
    open (my $reads_ofh, ">$reads_file") or die $!;
    

    my $read_counter = 0;

    ## get transcript structures:  if genome, can be complex. If transcriptome, should be simple (single coordinates per transcript).
    my $gene_obj_indexer_href = {};
    my $contig_to_gene_list_href;
    
    ## associate gene identifiers with contig id's.
    if ($gff3_file) {
        $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);
    }
    else {
        $contig_to_gene_list_href = &GTF_utils::index_GTF_gene_objs($gtf_file, $gene_obj_indexer_href);
    }
    

    foreach my $contig (keys %$contig_to_gene_list_href) {
        
        my $contig_seq = $transcript_seqs{$contig} or die "Error, no contig sequence for $contig";
        
        my @gene_ids = @{$contig_to_gene_list_href->{$contig}};
        
        foreach my $gene_id (@gene_ids) {
            
            print STDERR "\n-processing $contig :: $gene_id\n";
            
            my $gene_obj = $gene_obj_indexer_href->{$gene_id};

            
            						
			my $orientation = $gene_obj->get_orientation();
            

            my %iso_coordsets = &get_isoform_to_coordsets($gene_obj);
            my %iso_cdna_seqs = &get_isoform_seqs(\$contig_seq, \%iso_coordsets, $orientation);
            
            
            
            my @all_oriented_reads;

            foreach my $isoform_acc (keys %iso_coordsets) {
                
                my $coordset = $iso_coordsets{$isoform_acc} or die "Error, no coordset set for $isoform_acc";
                
                my $isoform_seq = $iso_cdna_seqs{$isoform_acc} or die "Error, no cdnaseq for $isoform_acc";
                print $cdna_ofh ">$isoform_acc\n$isoform_seq\n";
                
                my $read_simulator = new Simulate::Uniform_Read_Generator( { coordsets => $coordset,
                                                                             mean_fragment_length => $mean_frag_length,
                                                                             fragment_length_stdev => $frag_length_stdev,
                                                                             read_length => $read_length,
                                                                         }
                    );
                
                print STDERR "-simulating reads for $isoform_acc\n";
                
                my @paired_reads = $read_simulator->simulate_paired_reads_uniformly_across_seq();
                
                
                foreach my $paired_read (@paired_reads) {
                        
                    my ($left_read, $right_read) = @$paired_read;
                    
                    
                    if ($orientation eq '-') {
                        
                        ($left_read, $right_read) = ($right_read, $left_read);
                    }
                    

                    my $read_name = "R" . sprintf("%09s", ++$read_counter);
                    
                    
                    my @oriented_reads;
                    if ($SS_lib_type) {

                        #  (+) and (-) are respective to the genome, not to the transcript. 
                        
                        if ($SS_lib_type eq "RF") {
                            push (@oriented_reads, [$right_read, 1, '-', $read_name], [$left_read, 2, '+', $read_name]);
                        }
                        elsif ($SS_lib_type eq "FR") {
                            push (@oriented_reads, [$left_read, 1, '+', $read_name], [$right_read, 2, '-', $read_name]);
                        }
                        elsif ($SS_lib_type eq "F") {
                            push (@oriented_reads, [$left_read, 0, '+', $read_name]);
                        }
                        elsif ($SS_lib_type eq "R") {
                            push (@oriented_reads, [$right_read, 0, '-', $read_name]);
                        }
                        else {
                            die "Error, SS_lib_type [$SS_lib_type] not recognized";
                        }
                    }
                    else {
                        # not strand-specific
                        my ($pos_1, $pos_2) = shuffle(1,2);

                        if ($paired_flag) {
                            push (@oriented_reads, [$left_read, $pos_1, '+', $read_name], [$right_read, $pos_2, '-', $read_name]);

                        }
                        else {
                            # unpaired, choose one read or the other randomly.
                            my @read_options = shuffle ( [$left_read, 0, '+', $read_name], [$right_read, 0, '-', $read_name]);
                            push (@oriented_reads, $read_options[0]);
                        }
                    }
                    
                    my %opposite =  ( '+' => '-',
                                      '-' => '+',
                                      );

                    ## write reads to file
                    foreach my $oriented_read (@oriented_reads) {
                        my ($read, $frag_pos, $genome_orient, $name) = @$oriented_read;
                        my ($read_orient) = ($orientation eq '-') ? $opposite{$genome_orient} : $genome_orient;
                        
                        my $read_seq = &construct_seq_from_coordset($read, \$contig_seq, $read_orient);
                        my $fpos = ($frag_pos == 0) ? 1 : $frag_pos;
                        print $reads_ofh ">$name/$fpos\n$read_seq\n";
                    }
                    

                    push (@all_oriented_reads, [@oriented_reads]); # group pairs if paired.
                    
                } # end of foreach paired read
            } # end of foreach isoform
             

            
            ## examine isoform compatibility
            my %compatibility_map;
            
            foreach my $oriented_readset (@all_oriented_reads) {
                
                foreach my $isoform_acc (keys %iso_coordsets) {
                    
                    my $isoform_coordset = $iso_coordsets{$isoform_acc};
                    
                    if (  
                        ($paired_flag
                         && &Overlap_info::compatible_overlap_A_contains_B($isoform_coordset, $oriented_readset->[0]->[0])
                         && &Overlap_info::compatible_overlap_A_contains_B($isoform_coordset, $oriented_readset->[1]->[0])   # both reads must be compatible with isoform
                        )
                        ||
                        # unpaired
                        ( (! $paired_flag) && &Overlap_info::compatible_overlap_A_contains_B($isoform_coordset, $oriented_readset->[0]->[0]) ) ) {
                        
                        $compatibility_map{$oriented_readset}->{$isoform_acc} = 1;
                    }
                }
            } 
            # end of isoform compatibility map
            
                        
            # write reads to sam file
            my $read_counter = 0;
            foreach my $oriented_readset (@all_oriented_reads) {
                
                $read_counter++;
                print STDERR "\r[$read_counter] reads written. ";
                
                
                &write_readset_to_SAM_file($oriented_readset, $genome_sam_FH, $orientation, $contig, \$contig_seq);
                
                foreach my $isoform (keys %{$compatibility_map{$oriented_readset}}) {
                    
                    my $iso_coordset = $iso_coordsets{$isoform};
                    
                    my @isoform_oriented_readset;
                    
                    foreach my $read_info_aref (@$oriented_readset) {
                        
                        my @read_info = @$read_info_aref;
                        my $read_coordset = $read_info[0];
                        
                        my $isoform_adjusted_read_coordset = &transform_to_isoform_coordinates($read_coordset, $iso_coordset, $orientation);
                        $read_info[0] = $isoform_adjusted_read_coordset;
                        
                        push (@isoform_oriented_readset, [@read_info]);
                    }
                    
                    my $isoform_seq = $iso_cdna_seqs{$isoform} or die "Error, no isoform sequence for $isoform";;
                    
                    &write_readset_to_SAM_file(\@isoform_oriented_readset, $trans_sam_FH, '+', $isoform, \$isoform_seq);
                    
                    
                    
                } # end of writing isoform-level reads
                
                
                
            } # end of sam writing for genome reads

        } # end of each gene processing
        
    } # end of scaffold processing
    
    
    close $genome_sam_FH;
    close $trans_sam_FH;
    close $cdna_ofh;
    
    ###################################
    ## Process transcriptome SAM files
    ## write to sam and bam formats.
    
    ## prepare transcriptome bams
    my $cdna_fai = "$cdna_file.fai";
    my $cmd = "samtools faidx $cdna_file";
    &process_cmd($cmd);
    
    
#    $cmd = "samtools view -bt $cdna_fai $trans_sam_outfile | samtools sort -n - $trans_sam_outfile.nameSorted";
#    &process_cmd($cmd);
    
    $cmd = "samtools view -bt $cdna_fai $trans_sam_outfile | samtools sort - $trans_sam_outfile.coordSorted";
    &process_cmd($cmd);
    unlink("$trans_sam_outfile");
    
=strand_sep_trans

    if ($SS_lib_type) {
        $cmd = "$FindBin::RealBin/../support_scripts/SAM_strand_separator.pl $trans_sam_outfile.coordSorted.bam $SS_lib_type";
        &process_cmd($cmd);

        foreach my $sam_file ("$trans_sam_outfile.coordSorted.bam.+.sam", "$trans_sam_outfile.coordSorted.bam.-.sam") {

            if (-s $sam_file) {
                # convert to bam
                my $bam_file = $sam_file;
                $bam_file =~ s/sam$/bam/;
                
                my $cmd = "samtools view -bt $cdna_fai $sam_file > $bam_file";
                &process_cmd($cmd);
                unlink($sam_file);

                $cmd = "samtools index $bam_file";
                &process_cmd($cmd);
            }
        }
    }

=cut

    ##############################
    ## Process genome file
    ## sort genome file
            
    # prepare genome bams
    my $genome_fai = "$genome_file.fai";
    unless (-s $genome_fai) {
        my $cmd = "samtools faidx $genome_file";
        &process_cmd($cmd);
    }
    
    $cmd = "samtools view -bt $genome_fai $genome_sam_outfile | samtools sort - $genome_sam_outfile.coordSorted";
    &process_cmd($cmd);
    
    unlink("$genome_sam_outfile");

    $cmd = "samtools index $genome_sam_outfile.coordSorted.bam";
    &process_cmd($cmd);

    
=strand_sep_genome

    if ($SS_lib_type) {
        $cmd = "$FindBin::RealBin/../support_scripts/SAM_strand_separator.pl $genome_sam_outfile.coordSorted.bam $SS_lib_type";
        &process_cmd($cmd);

        foreach my $sam_file ("$genome_sam_outfile.coordSorted.bam.+.sam", "$genome_sam_outfile.coordSorted.bam.-.sam") {

            if (-s $sam_file) {
                
                # convert to bam
                my $bam_file = $sam_file;
                $bam_file =~ s/sam$/bam/;
                
                my $cmd = "samtools view -bt $genome_fai $sam_file > $bam_file";
                &process_cmd($cmd);
                
                $cmd = "samtools index $bam_file";
                &process_cmd($cmd);
            }
            
            unlink("$sam_file");
            
        }
    }
    
=cut

    exit(0);
}



####
sub compute_transcript_length {
    my ($coordset) = @_;

    my $len = 0;
    foreach my $coord_pair (@$coordset) {
        my ($lend, $rend) = @$coord_pair;
        $len += $rend - $lend + 1;
    }

    return($len);
}



####
sub get_isoform_to_coordsets {
    my ($gene_obj) = @_;
    
    my %isoform_to_coordsets;
    
    foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {
        
        my $isoform_id = $isoform->{Model_feat_name};
        
        my @exons = $isoform->get_exons();
        
        my @coords;
        
        foreach my $exon (@exons) {
            
            my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
            
            push (@coords, [$lend, $rend]);
        }
        
        @coords = sort {$a->[0]<=>$b->[0]} @coords;
        
        $isoform_to_coordsets{$isoform_id} = \@coords;
    }
    
    
    return(%isoform_to_coordsets);
    
}


####
sub construct_seq_from_coordset {
    my ($read, $genome_sref, $orientation) = @_;
    
    my $read_seq = "";
    foreach my $coords (@$read) {

        my ($lend, $rend) = @$coords;
    
        if ($rend > length($$genome_sref)) {
            confess "Error, $rend > " . length($$genome_sref) . ", " . Dumper($read);
        }
    
        my $seq = substr($$genome_sref, $lend-1, $rend - $lend + 1);
        $read_seq .= $seq;
    }

	if ($orientation eq '-') {
		$read_seq = &reverse_complement($read_seq);
	}
	
    return($read_seq);
}


####
sub get_cdna_coords {
    my ($read) = @_;

    my $last_cdna_pos = 0;

    my @cdna_coords;
    
    foreach my $coordset (@$read) {
        my ($lend, $rend) = @$coordset;
        my $left_cdna = $last_cdna_pos + 1;
        my $right_cdna = $last_cdna_pos + $rend - $lend + 1;
        push (@cdna_coords, [$left_cdna, $right_cdna]);
        
        $last_cdna_pos = $right_cdna;
    }

    return(@cdna_coords);
}


####
sub create_SAM_entry {
    my ($read_name, $read_coordset, $frag_pos, $target_name, $read_seq, $orientation) = @_;
	
    my @fields;
    $#fields = 10;
    
    $fields[0] = $read_name;
    $fields[1] = 0;
    $fields[2] = $target_name;
    $fields[3] = $read_coordset->[0]->[0]; # first coordinate
    $fields[4] = 255;
    $fields[9] = $read_seq;
    $fields[10] = 'B' x length($read_seq);
    
    my @cdna_coords = &get_cdna_coords($read_coordset);
    
    my $cigar_string = &CIGAR::construct_cigar($read_coordset, \@cdna_coords, length($read_seq));
	$cigar_string =~ s/D/N/g;
	
    $fields[5] = $cigar_string;
    $fields[6] = '*';
	$fields[7] = 0;
	$fields[8] = 0;
    
    my $sam = new SAM_entry(join("\t", @fields));
    
    $sam->set_query_strand($orientation);

    if ($frag_pos == 1) {
        $sam->set_first_in_pair(1);
    }
    elsif ($frag_pos == 2) {
        $sam->set_second_in_pair(2);
    }
    
    
    return ($sam);
}



####
sub write_SAM_transcriptome {
    my ($read_name, $read_coordset, $target_name, $FH, $read_seq, $orientation, $cdna_seq, $isoform_coordset) = @_;
	
    my @fields;
    $#fields = 10;
    
    $fields[0] = $read_name;
    $fields[1] = 0;
    $fields[2] = $target_name;

	if ($orientation eq '-') {
		$read_seq = &reverse_complement($read_seq);
	}
	
	my $start_pos = index($cdna_seq, $read_seq);
	if ($start_pos < 0) {
		print STDERR "Read: " . Dumper($read_coordset) . "Transcript: " . Dumper($isoform_coordset);
		confess "Error, cannot map read sequence to cdna:\nRead: $read_seq\nCDNA: $cdna_seq\n";
	}
	
	$fields[3] = $start_pos + 1;
	
    $fields[4] = 255;
    $fields[9] = $read_seq;
    $fields[10] = 'B' x length($read_seq);
    
    my @cdna_coords = &get_cdna_coords($read_coordset);
    
    my $cigar_string = length($read_seq) . "M";
    $fields[5] = $cigar_string;
    $fields[6] = '*';
	$fields[7] = 0;
	$fields[8] = 0;
    
    my $sam = new SAM_entry(join("\t", @fields));
    my $sam_line = $sam->toString();
	$sam_line .= "\tXS:A:+";

    print $FH $sam_line . "\n";



    return;
}

####
sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";
	my $ret = system($cmd);
	if ($ret) {
		confess "Error, CMD: $cmd died with ret $ret";
	}

	return;
}


####
sub parse_frag_counts {
    my ($file) = @_;

    my %counts;

    open (my $fh, $file) or die "Error, cannot open file $file";
    while(<$fh>) {
        chomp;
        my ($acc, $count) = split(/\t/);
        unless ($count =~ /^\d+$/) {
            die "Error, cannot parse count from $file, entry: $_";
        }
        $counts{$acc} = $count;
    }

    return(%counts);
}




####
sub write_readset_to_SAM_file {
    my ($oriented_readset, $fh, $transcribed_orientation, $contig, $contig_seq_sref) = @_;
        
    my @sam_entries;
    
    foreach my $oriented_read (@$oriented_readset) {
        
        my ($read, $frag_pos, $read_orient, $read_name) = @$oriented_read;
        my $read_seq = &construct_seq_from_coordset($read, $contig_seq_sref, '+'); # reads always in plus orientation for genome.
        
        if ($transcribed_orientation eq '-') {
            ## flip the read orientation
            $read_orient = ($read_orient eq '+') ? '-' : '+';
        }
        
        my $sam_entry = &create_SAM_entry($read_name, $read, $frag_pos, $contig, $read_seq, $read_orient);
        push (@sam_entries, $sam_entry);
        
    }
    
    if ($paired_flag) {
        $sam_entries[0]->set_paired(1);
        $sam_entries[1]->set_paired(1);
        
        $sam_entries[0]->set_proper_pair(1);
        $sam_entries[0]->set_mate_strand( $sam_entries[1]->get_query_strand());
        $sam_entries[0]->set_mate_scaffold_name( $sam_entries[1]->get_scaffold_name());
        $sam_entries[0]->set_mate_scaffold_position( $sam_entries[1]->get_scaffold_position());
        
        $sam_entries[1]->set_proper_pair(1);
        $sam_entries[1]->set_mate_strand( $sam_entries[0]->get_query_strand());
        $sam_entries[1]->set_mate_scaffold_name( $sam_entries[0]->get_scaffold_name());
        $sam_entries[1]->set_mate_scaffold_position( $sam_entries[0]->get_scaffold_position());
    }
    
    foreach my $sam_entry (@sam_entries) {
        print $fh $sam_entry->toString();
        if ($SS_lib_type) {
            my $XS_flag = "XS:A:$transcribed_orientation";
            print $fh "\t$XS_flag";
        }
        
        print $fh "\n";
    }

    return;
}


####
sub transform_to_isoform_coordinates {
    my ($read_coordset, $iso_coordset, $orientation) = @_;
    
    my @iso_coord_structs;

    @$iso_coordset = sort {$a->[0]<=>$b->[0]} @$iso_coordset;  # sort by left segment boundary coordinate

    
    my $cdna_len = 0;
    
    foreach my $iso_segment (@$iso_coordset) {
        my ($lend, $rend) = @$iso_segment;

        my $cdna_lend = $cdna_len + 1;
        my $cdna_rend = $cdna_len + $rend - $lend + 1;
        
        $cdna_len += $rend - $lend + 1;

        push (@iso_coord_structs, { genome_lend => $lend,
                                    genome_rend => $rend,
                                    
                                    cdna_lend => $cdna_lend,
                                    cdna_rend => $cdna_rend,
                                });

    }

    my ($genome_left_bound, $genome_right_bound) = ($read_coordset->[0]->[0], $read_coordset->[$#$read_coordset]->[1]);

    my $cdna_left_bound = &_convert_genome_to_cdna_coord($genome_left_bound, \@iso_coord_structs);
    my $cdna_right_bound = &_convert_genome_to_cdna_coord($genome_right_bound, \@iso_coord_structs);

    if ($orientation eq '-') {
        # revcomp the coordinates

        ($cdna_left_bound, $cdna_right_bound) = ($cdna_len - $cdna_right_bound + 1, $cdna_len - $cdna_left_bound + 1);
    }

    return ([ [$cdna_left_bound, $cdna_right_bound] ]);  # single segment
}

####
sub _convert_genome_to_cdna_coord {
    my ($genome_coord, $iso_coord_structs) = @_;


    foreach my $struct (@$iso_coord_structs) {
        
        if ($genome_coord >= $struct->{genome_lend} && $genome_coord <= $struct->{genome_rend}) {

            my $cdna_coord = $struct->{cdna_lend} + $genome_coord - $struct->{genome_lend};

            return($cdna_coord);
        }
    }

    die "Error, couldn't map genome coordinate ($genome_coord) to iso coordinate within: " . Dumper($iso_coord_structs);

}


####
sub get_isoform_seqs {
    my ($genome_seq_sref, $iso_coordsets_href, $orientation) = @_;

    my %cdna_seqs;

    foreach my $isoform_acc (keys %$iso_coordsets_href) {

        my $iso_coordset = $iso_coordsets_href->{$isoform_acc};

        my $cdna_seq = &construct_seq_from_coordset($iso_coordset, $genome_seq_sref, $orientation); 
    
        $cdna_seqs{$isoform_acc} = $cdna_seq;
    }

    return(%cdna_seqs);
}
