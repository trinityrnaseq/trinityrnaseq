#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;


my $usage = "usage: $0 genes.gff3 genome.fasta\n\n";

my $genes_gff3 = $ARGV[0] or die $usage;
my $genome_fasta_file = $ARGV[1] or die $usage;


my $fasta_reader = new Fasta_reader($genome_fasta_file);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($genes_gff3, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = $genome{$asmbl_id} or die "Error, cannot find sequence for $asmbl_id"; #cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
        my $orientation = $gene_obj_ref->get_orientation();
        
        my $intron_text = "";
        
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
 
            my @intron_coords = $isoform->get_intron_coordinates();
            
            
            my $trans_id = $isoform->{Model_feat_name};
            
            foreach my $intron (@intron_coords) {
                my ($end5, $end3) = @$intron;
                my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
                my $left_splice_dinuc = substr($genome_seq, $lend-1, 2);
                my $right_splice_dinuc = substr($genome_seq, $rend-1-1, 2);
                
                $intron_text .= join("\t", $gene_id, $trans_id, $asmbl_id, "$lend-$rend", $orientation, "$left_splice_dinuc..$right_splice_dinuc") . "\n";
            }
        }
        
        if ($intron_text) {
            print "$intron_text\n";
        }
        
        
            
        
    }
}


exit(0);

