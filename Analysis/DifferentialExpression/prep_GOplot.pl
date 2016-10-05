#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use DelimParser;

my $usage = "\n\n\tusage: $0 Trinotate_report.xls.gene_ontology DE_analysis.GOseq.enriched\n\n";

my $gene_ontology_assignments_file = $ARGV[0] or die $usage;
my $enriched_file = $ARGV[1] or die $usage;


main: {

     my %enriched_GO = &parse_GO_enriched($enriched_file);
     
     my %genes_with_enriched_GO = &parse_GO_assignments($gene_ontology_assignments_file, \%enriched_GO);

     my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", ['Category', 'ID', 'Term', 'Genes', 'adj_pval']);

     
     foreach my $row (values %enriched_GO) {

         # format wanted:
         #  Category         ID              Term
         #1       BP GO:0007507 heart development
         #    Genes
         #    1 DLC1, NRP2, NRP1, EDN1, PDLIM3, GJA1, TTN, GJA5, ZIC3, TGFB2, CERKL, GATA6, COL4A3BP, GAB1, SEMA3C, MKL2, SLC22A5, MB, PTPRJ, RXRA, VANGL2, MYH6, TNNT2, HHEX, MURC, MIB1, FOXC2, FOXC1, ADAM19, MYL2, TCAP, EGLN1, SOX9, ITGB1, CHD7, HEXIM1, PKD2, NFATC4, PCSK5, ACTC1, TGFBR2, NF1, HSPG2, SMAD3, TBX1, TNNI3, CSRP3, FOXP1, KCNJ8, PLN, TSC2, ATP6V0A1, TGFBR3, HDAC9
         #    adj_pval
         #    1 2.17e-06

         my $go_term_info = $row->{go_term};
         my ($go_type, $go_descr) = split(/\s+/, $go_term_info, 2);
         my $go_id = $row->{'category'};
         
         $tab_writer->write_row( { 'Category' => $go_type,
                                   'ID' => $go_id,
                                   'Term' => $go_descr,
                                   'Genes' => $genes_with_enriched_GO{$go_id},
                                   'adj_pval' => $row->{over_represented_FDR},
                                 } );
         
     }
     

     exit(0);
}

####
sub parse_GO_enriched {
    my ($enriched_file) = @_;

    open (my $fh, $enriched_file) or die "Error, cannot open file $enriched_file";
    my $tab_reader = new DelimParser::Reader($fh, "\t");

    my %enriched_GO;
    
    while (my $row = $tab_reader->get_row()) {
        
        my $category = $row->{category} or die "Error, cannot identify 'category' value";
        $enriched_GO{$category} = $row;
    }

    return(%enriched_GO);
}

####
sub parse_GO_assignments {
    my ($gene_ontology_assignment_file, $enriched_GO_href) = @_;

    my %GO_to_gene_ids;
    
    open (my $fh, $gene_ontology_assignment_file) or die "Error, cannot open file $gene_ontology_assignment_file";
    while (<$fh>) {
        chomp;
        my ($gene_id, $GO_assignments) = split(/\t/);
        my @GO = split(/,/, $GO_assignments);

        foreach my $GO_id (@GO) {
            if (exists $enriched_GO_href->{$GO_id}) {
                $GO_to_gene_ids{$GO_id}->{$gene_id} = 1;
            }
        }
    }
    close $fh;

    my %genes_enriched_GO;

    foreach my $GO_id (keys %GO_to_gene_ids) {
        my $gene_ids_href = $GO_to_gene_ids{$GO_id};
        my @gene_list = sort keys %$gene_ids_href;
        $genes_enriched_GO{$GO_id} = join(", ", @gene_list);
    }

    return(%genes_enriched_GO);

}


        
