#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$ENV{TRINITY_HOME}/PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;

my $usage = "\n\nusage: $0 gff3_file genome_db [flank]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $flank = $ARGV[2] || 0;

my $fasta_reader = new Fasta_reader($fasta_db);
my %genome = $fasta_reader->retrieve_all_seqs_hash();

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

my $outdir = "gene_contigs";
mkdir $outdir or die "Error, $outdir already exists";
my $gene_contigs_per_bin = 100;
my $gene_counter = 0;


foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my $genome_seq = $genome{$asmbl_id} or die "Error, cannot find sequence for $asmbl_id"; #cdbyank_linear($asmbl_id, $fasta_db);
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {

        my $gene_id_for_filename = $gene_id;
        $gene_id_for_filename =~ s/\W/_/g;
        
        $gene_counter++;
        my $bindir = "$outdir/bin_" . int($gene_counter/$gene_contigs_per_bin) . "/$gene_id_for_filename";
        print STDERR "[$gene_counter] -processing $bindir\n";
        &process_cmd("mkdir -p $bindir");
        
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		
        my ($gene_contig_lend, $gene_contig_rend) = sort {$a<=>$b} $gene_obj_ref->get_gene_span();
        
        $gene_contig_lend -= $flank;
        $gene_contig_rend += $flank;

        my $contig_filename = "$bindir/gene.fa";
        open (my $ofh, ">$contig_filename") or die $!;
        my $gene_seq = substr($genome_seq, $gene_contig_lend-1, $gene_contig_rend - $gene_contig_lend + 1);
        print $ofh ">$gene_id\n$gene_seq\n";
        close $ofh;

        $gene_obj_ref->adjust_gene_coordinates(-1 * ($gene_contig_lend-1));
        
        # update gene identifiers.
        foreach my $iso_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
            $iso_obj->{asmbl_id} = $gene_id;
        }
        

        my $gff3_filename = "$bindir/gene.gff3";
        open ($ofh, ">$gff3_filename") or die $!;
        print $ofh $gene_obj_ref->to_GFF3_format();
        close $ofh;
        

        #if ($gene_counter > 10) { last; }
        
    }
}


exit(0);

####
sub process_cmd {
    my ($cmd) = @_;

    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }
    
    return;
}
