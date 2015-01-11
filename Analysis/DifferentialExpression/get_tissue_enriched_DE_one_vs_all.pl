#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 MAX_FDR MIN_FC\n\n"
    . " \n\tRun this in the edgeR results directory.\n\n";

my $MAX_FDR = $ARGV[0] or die $usage;
my $MIN_FC = $ARGV[1] or die $usage;

my $MIN_LOGFC = log($MIN_FC)/log(2);

main: {

    my @DE_files = <*.DE_results>;
    
    my %sample_to_comparison_count;
    my %sample_to_DE_genes;
    
    foreach my $DE_file (@DE_files) {
        
        $DE_file =~ /\.([^\.]+)_vs_([^\.]+)/ or die "Error, cannot extract sample names from filename: $DE_file";
        
        my $sample_A = $1;
        my $sample_B = $2;
        
        $sample_to_comparison_count{$sample_A}++;
        $sample_to_comparison_count{$sample_B}++;
        
        open (my $fh, $DE_file) or die "Error, cannot open file $DE_file";
        my $header = <$fh>;
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $logFC = $x[1];
            my $FDR = $x[4];
            
            unless ($FDR <= $MAX_FDR && abs($logFC) >= $MIN_LOGFC) { next; }
            
            my $gene = $x[0];
            
            if ($logFC < 0) {
                $sample_to_DE_genes{$sample_A}->{$gene}++;
            }
            else {
                $sample_to_DE_genes{$sample_B}->{$gene}++;
            }
        }
        close $fh;
        
    }
    
    foreach my $sample (keys %sample_to_DE_genes) {
        foreach my $gene (keys %{$sample_to_DE_genes{$sample}}) {

            my $count_DE_induced = $sample_to_DE_genes{$sample}->{$gene};
            
            my $num_sample_comparisons = $sample_to_comparison_count{$sample};
            if ($count_DE_induced == $num_sample_comparisons) {
                ## DE and induced in each comparison
                
                print "$gene\t$sample\n";
            }
        }
    }

    exit(0);
}
