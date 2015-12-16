#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $help_flag;

my $usage = <<__EOUSAGE__;


####################################################################
#
#  Required:
#
#  --GO_annots <string>            GO assignments to genes.
#
#  --gene_factors_file <string>    gene factors file
#
#  --factor <string>               the gene factor defining the genes of interest.
#
#  --GO_id <string>                GO ID to examine.
#
#  --fpkm_matrix <string>          FPKM matrix file.
#
#  Optional:
#
#  --output_prefix <string>        output prefix  (default: 'factor,GO_id')
#
#  --samples <string>              samples description file (for heatmaps)
#
####################################################################

__EOUSAGE__

    ;

my $GO_annots_file;
my $gene_factors_file;
my $factor;
my $GO_id;
my $fpkm_matrix;
my $output_prefix;
my $samples_file;

&GetOptions ( 'h' => \$help_flag,
              'GO_annots=s' => \$GO_annots_file,
              'gene_factors_file=s' => \$gene_factors_file,
              'factor=s' => \$factor,
              'GO_id=s' => \$GO_id,
              'fpkm_matrix=s' => \$fpkm_matrix,
              'output_prefix=s' => \$output_prefix,
              'samples=s' => \$samples_file,
    );


if ($help_flag) { 
    die $usage;
}

unless ($GO_annots_file && $gene_factors_file && $factor && $GO_id && $fpkm_matrix) {
    die $usage;
}

unless ($output_prefix) {
    $output_prefix = join(",", $factor, $GO_id);
}

main: {

    my %genes = &get_genes_assigned_to_factor($gene_factors_file, $factor);
    
    my %genes_with_GO = &get_subset_of_genes_with_GO_id($GO_annots_file, $GO_id, \%genes);

    my $fpkm_outfile = "$output_prefix.fpkm";
    open (my $ofh, ">$fpkm_outfile") or die "Error, cannot write to $fpkm_outfile";
    open (my $fh, $fpkm_matrix) or die "Error, cannot open file $fpkm_matrix";
    my $header = <$fh>;
    print $ofh $header;
    
    my %genes_want = %genes_with_GO;
    while (<$fh>) {
        my $line = $_;
        my @x = split(/\t/);
        my $acc = $x[0];
        if ($genes_with_GO{$acc}) {
            print $ofh $line;
            delete $genes_want{$acc};
        }
    }
    close $ofh;

    if (%genes_want) {
        use Data::Dumper;
        die "Error, didn't find fpkm matrix rows for genes: " . Dumper(\%genes_want);
    }
    
    ## generate a heatmap
    my $cmd = "$FindBin::RealBin/PtR -m $fpkm_outfile --log2 --heatmap --gene_dist euclidean --sample_dist euclidean";
    if ($samples_file) {
        $cmd .= " -s $samples_file ";
    }
    
    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }
    
    exit(0);
    
}

####
sub get_genes_assigned_to_factor {
    my ($gene_factors_file, $factor_want) = @_;

    my %genes;

    open (my $fh, $gene_factors_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($gene, $factor) = split(/\t/);
        if ($factor eq $factor_want) {
            
            $genes{$gene} = 1;
        }
    }
    close $fh;

    unless (%genes) {
        confess "Error, no genes found assigned to factor [$factor_want] in file $gene_factors_file";
    }

    return(%genes);
}


####
sub get_subset_of_genes_with_GO_id {
    my ($GO_annots_file, $GO_id, $genes_href) = @_;

    my %genes_with_GO;

    open (my $fh, $GO_annots_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($gene, $go_annots) = split(/\t/);
        if ($genes_href->{$gene} && $go_annots =~ /$GO_id/) {
            $genes_with_GO{$gene} = 1;
        }
    }
    close $fh;

    unless (%genes_with_GO) {
        confess "Error, no genes found with GO: [$GO_id] in file: $GO_annots_file for genes: " . Dumper($genes_href);
    }

    return(%genes_with_GO);
}


