#!/usr/bin/env perl

package GFF3_alignment_utils;

use strict;
use warnings;
use Carp;

use Gene_obj;
use Gene_obj_indexer;
use CDNA::Alignment_segment;
use CDNA::CDNA_alignment;
use File::Basename;

__run_test() unless caller;

sub index_alignment_objs {
    my ($gff3_alignment_file, $genome_alignment_indexer_href) = @_;
    
    unless ($gff3_alignment_file && -s $gff3_alignment_file) {
        confess "Error, cannot find or open file $gff3_alignment_file";
    }
    unless (ref $genome_alignment_indexer_href) {
        confess "Error, need genome indexer href as param ";
    }
    
       	
	my %genome_trans_to_alignment_segments;
    my %trans_to_gene_id;

	
	open (my $fh, $gff3_alignment_file) or die "Error, cannot open file $gff3_alignment_file";
	while (<$fh>) {
		chomp;
		
		unless (/\w/) { next; }
		
		my @x = split(/\t/);

		unless (scalar (@x) >= 8 && $x[8] =~ /ID=/) {
			print STDERR "ignoring line: $_\n";
			next;
		}
		
		my $scaff = $x[0];
		my $type = $x[2];
		my $lend = $x[3];
		my $rend = $x[4];
        my $per_id = $x[5];
        if ($per_id eq ".") { $per_id = 100; } # making an assumption here.
        
		my $orient = $x[6];
		
		my $info = $x[8];
		
		my @parts = split(/;/, $info);
		my %atts;
		foreach my $part (@parts) {
			$part =~ s/^\s+|\s+$//;
			$part =~ s/\"//g;
			my ($att, $val) = split(/=/, $part);
			
			if (exists $atts{$att}) {
				die "Error, already defined attribute $att in $_";
			}
			
			$atts{$att} = $val;
		}

		my $gene_id = $atts{ID} or die "Error, no gene_id at $_";
		my $trans_id = $atts{Target} or die "Error, no trans_id at $_";
		{
			my @pieces = split(/\s+/, $trans_id);
			$trans_id = shift @pieces;
		}
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        $info =~ /Target=\S+ (\d+) (\d+)/ or die "Error, cannot extract match coordinates from info: $info";
        my $cdna_seg_lend = $1;
        my $cdna_seg_rend = $2;
        
        ($cdna_seg_lend, $cdna_seg_rend) = sort {$a<=>$b} ($cdna_seg_lend, $cdna_seg_rend); # always + orient for transcript coords.
        
        
        my $alignment_segment = new CDNA::Alignment_segment($end5, $end3, $cdna_seg_lend, $cdna_seg_rend, $per_id);
        

        push (@{$genome_trans_to_alignment_segments{$scaff}->{$trans_id}}, $alignment_segment);
        
		$trans_to_gene_id{$trans_id} = $gene_id;
        
	}
    

    my %scaff_to_align_list;
    
    
	## Output genes in gff3 format:

	foreach my $scaff (sort keys %genome_trans_to_alignment_segments) {
        
        my @alignment_accs = keys %{$genome_trans_to_alignment_segments{$scaff}};

        foreach my $alignment_acc (@alignment_accs) {

            my $segments_aref = $genome_trans_to_alignment_segments{$scaff}->{$alignment_acc};
            
            ## determine cdna length
            my @cdna_coords;
            foreach my $segment (@$segments_aref) {
                push (@cdna_coords, $segment->get_mcoords());
            }
            @cdna_coords = sort {$a<=>$b} @cdna_coords;
            my $max_coord = pop @cdna_coords;

            my $cdna_alignment_obj = new CDNA::CDNA_alignment($max_coord, $segments_aref);
            $cdna_alignment_obj->set_acc($alignment_acc);
            $cdna_alignment_obj->{genome_acc} = $scaff;
            
            my $gene_id = $trans_to_gene_id{$alignment_acc} or confess "Error no gene_id for acc: $alignment_acc";

            $cdna_alignment_obj->{gene_id} = $gene_id;
            

            $cdna_alignment_obj->{source} = basename($gff3_alignment_file);
            
            if (ref $genome_alignment_indexer_href eq "Gene_obj_indexer") {
                $genome_alignment_indexer_href->store_gene($alignment_acc, $cdna_alignment_obj);
                
            }
            else {
                                
                $genome_alignment_indexer_href->{$alignment_acc} = $cdna_alignment_obj;
            }
            push (@{$scaff_to_align_list{$scaff}}, $alignment_acc);
        }
    }

    return(%scaff_to_align_list);
}



#################
## Testing
#################


sub __run_test {
    
    my $usage = "usage: $0 file.alignment.gff3\n\n";

    my $gff3_file = $ARGV[0] or die $usage;

    my $indexer = {};
    my %scaff_to_alignments = &index_alignment_objs($gff3_file, $indexer);
    
    
    foreach my $scaffold (keys %scaff_to_alignments) {

        my @align_ids = @{$scaff_to_alignments{$scaffold}};

        foreach my $align_id (@align_ids) {
            my $cdna_obj = $indexer->{$align_id};

            print $cdna_obj->toString();
        }
    }
    
    
    
    exit(0);
    


}







1; #EOM

