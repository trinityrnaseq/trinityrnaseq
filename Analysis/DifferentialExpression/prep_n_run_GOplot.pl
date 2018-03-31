#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use DelimParser;
use File::Basename;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $help_flag;

my $usage = <<__EOUSAGE__;

###################################################################
#
# --GO_annots <string>         Trinotate_report.xls.gene_ontology (genes or transcripts - be consistent w/ DE analysis being done here)
#
# --DE_subset <string>         The 'DE_analysis.DE_subset' file.
#
# --DE_GO_enriched <string>    The 'DE.subset.GOseq.enriched' file.
#
# --tmpdir <string>            Path to put prep files for GOplot
#
# --pdf_filename <string>      filename for pdf plot output file
#
###################################################################


__EOUSAGE__

    ;


my $gene_ontology_assignments_file;
my $enriched_file;
my $DE_results_file;
my $tmpdir;
my $pdf_filename;

&GetOptions ( 'h' => \$help_flag,
              'GO_annots=s' => \$gene_ontology_assignments_file,
              'DE_subset=s' => \$DE_results_file,
              'DE_GO_enriched=s' => \$enriched_file,
              'pdf_filename=s' => \$pdf_filename,
              'tmpdir=s' => \$tmpdir,
    );

if ($help_flag) {
    die $usage;
}

unless ($gene_ontology_assignments_file && $enriched_file && $DE_results_file && $tmpdir && $pdf_filename) {
    die $usage;
}



main: {

    if (! -d $tmpdir) {
        &process_cmd("mkdir -p $tmpdir");
    }
    
    
    my %enriched_GO = &parse_GO_enriched($enriched_file);
    
    my %DE_genes = &parse_DE_genes($DE_results_file);
    
    my %genes_with_enriched_GO = &parse_GO_assignments($gene_ontology_assignments_file, \%enriched_GO, \%DE_genes);



    ## Write the EC.david file:
    my $EC_david_file = "$tmpdir/EC.david";
    {
        open (my $ofh, ">$EC_david_file") or die "Error, cannot write to $EC_david_file";
        
        my $tab_writer = new DelimParser::Writer($ofh, "\t", ['Category', 'ID', 'Term', 'Genes', 'adj_pval']);
        
        foreach my $row (values %enriched_GO) {
            
            # format wanted:
            #  Category         ID              Term
            #1       BP GO:0007507 heart development
            #    Genes
            #    1 DLC1, NRP2, NRP1, EDN1, PDLIM3, GJA1, TTN, GJA5, ZIC3, TGFB2, CERKL, GATA6, COL4A3BP, GAB1, SEMA3C, MKL2, SLC22A5, MB, PTPRJ, RXRA, VANGL2, MYH6, TNNT2, HHEX, MURC, MIB1, FOXC2, FOXC1, ADAM19, MYL2, TCAP, EGLN1, SOX9, ITGB1, CHD7, HEXIM1, PKD2, NFATC4, PCSK5, ACTC1, TGFBR2, NF1, HSPG2, SMAD3, TBX1, TNNI3, CSRP3, FOXP1, KCNJ8, PLN, TSC2, ATP6V0A1, TGFBR3, HDAC9
            #    adj_pval
            #    1 2.17e-06
            
            my $go_term_info = $row->{go_term};
            
            my $go_id = $row->{'category'};

            if ($go_term_info eq "none") {
                print STDERR "WARNING, no GO term info found for: $go_id, skipping...\n";
                next;
            }
            
            my ($go_type, $go_descr) = split(/\s+/, $go_term_info, 2);
            
            my $genes = $genes_with_enriched_GO{$go_id};
            
            unless ($genes) {
                die "Error, no genes extracted for GO category: $go_id $go_type $go_descr ";
            }

            eval {
                $tab_writer->write_row( { 'Category' => $go_type,
                                          'ID' => $go_id,
                                          'Term' => $go_descr,
                                          'Genes' => $genes,
                                          'adj_pval' => $row->{over_represented_FDR},
                                        } );
            };
            if ($@) {
                confess "$@\n" . "row: " . Dumper($row);
            }
        }
        close $ofh;
    }

    ## write the EC.genelist file:
    my $EC_genelist_file = "$tmpdir/EC.genelist";
    {
        open (my $ofh, ">$EC_genelist_file") or die $!;

        my $tab_writer = new DelimParser::Writer($ofh, "\t", ['ID', 'logFC', 'adj.P.Val']);
        
        foreach my $DE_info_href (sort {$a->{'adj.P.Val'}<=>$b->{'adj.P.Val'}} values %DE_genes) {

            $tab_writer->write_row( { 'ID' => $DE_info_href->{ID},
                                      'logFC' => $DE_info_href->{logFC},
                                      'adj.P.Val' => $DE_info_href->{'adj.P.Val'},
                                    }
                );

        }

        close $ofh;
    }

    
    my $cmd = "$FindBin::RealBin/GOplot.Rscript --EC_david $EC_david_file --EC_genelist $EC_genelist_file --pdf_outfile $pdf_filename";
    &process_cmd($cmd);
    
    exit(0);
}


####
sub parse_DE_genes {
    my ($file) = @_;
    
    my %genes;
    
    open (my $fh, $file) or die "Error, cannot open file $file";

    my $header = <$fh>;
    unless ($header =~ /^sample/) {
        die "Error, file: $file has unexpected format... no 'sample' starting header.";
    }
    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $gene_id = $x[0];
        my $logFC = $x[3];
        my $fdr = $x[6];
        
        $genes{$gene_id} = { ID => $gene_id,
                             logFC => $logFC,
                             'adj.P.Val' => $fdr,
        };
        
    }
    
    return(%genes);
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
    my ($gene_ontology_assignment_file, $enriched_GO_href, $DE_genes_href) = @_;

    my %GO_to_gene_ids;
    
    open (my $fh, $gene_ontology_assignment_file) or die "Error, cannot open file $gene_ontology_assignment_file";
    while (<$fh>) {
        chomp;
        my ($gene_id, $GO_assignments) = split(/\t/);

        unless (exists $DE_genes_href->{$gene_id}) { next; }
        
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

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";

    my $ret = system($cmd);

    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}

        
