#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;

my $usage = "\n\nusage: $0 gff3_file genome_fasta [flank=0]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $flank = $ARGV[2] || 0;

my ($upstream_flank, $downstream_flank) = (0,0);

if ($flank) {
	if ($flank =~ /:/) {
		($upstream_flank, $downstream_flank) = split (/:/, $flank);
	}
	else {
		($upstream_flank, $downstream_flank) = ($flank, $flank);
	}
}

if ($upstream_flank < 0 || $downstream_flank < 0) {
	die $usage;
}



my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = $genome{$asmbl_id} or die "Error, cannot find sequence for $asmbl_id"; #cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
 
            my $orientation = $isoform->get_orientation();
			my ($model_lend, $model_rend) = sort {$a<=>$b} $isoform->get_model_span();
			my ($gene_lend, $gene_rend) = sort {$a<=>$b} $isoform->get_gene_span();
			
            my $isoform_id = $isoform->{Model_feat_name};
            
            my $seq = $isoform->create_cDNA_sequence(\$genome_seq);
            if ($upstream_flank || $downstream_flank) {
                $seq = &add_flank($seq, $upstream_flank, $downstream_flank, $gene_lend, $gene_rend, $orientation, \$genome_seq);
            }
			
            unless ($seq) {
                print STDERR "-warning, no cDNA sequence for $isoform_id\n";
                next;
            }
            
            $seq =~ s/(\S{60})/$1\n/g; # make fasta format
            chomp $seq;
            
            my $com_name = $isoform->{com_name} || "";
            
            print ">$gene_id" . "::" . "$isoform_id\n$seq\n";
        }
    }
}


exit(0);


####
sub add_flank {
	my ($seq, $upstream_flank, $downstream_flank, $lend, $rend, $orientation, $genome_seq_ref) = @_;
	
	my $far_left = ($orientation eq '+') ? $lend - $upstream_flank : $lend - $downstream_flank;
	
	if ($far_left < 1) { $far_left = 1; }
	
	my $flank_right = ($orientation eq '+') ? $downstream_flank : $upstream_flank;

	my $left_seq = substr($$genome_seq_ref, $far_left - 1, $lend - $far_left);

	my $right_seq = substr($$genome_seq_ref, $rend, $flank_right);
	
	if ($orientation eq '+') {
		return (lc($left_seq) . uc($seq) . lc($right_seq));
	}
	else {
		return (lc(&reverse_complement($right_seq)) . uc($seq) . lc(&reverse_complement($left_seq)));
	}
}


