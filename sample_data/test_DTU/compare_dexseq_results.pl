#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 genome.dexseq  superT.dexseq\n\n";

my $genome_dexseq = $ARGV[0] or die $usage;
my $superT_dexseq = $ARGV[1] or die $usage;


main: {

    my %genome_signif = &parse_signif_scores($genome_dexseq);
    my %superT_signif = &parse_signif_scores($superT_dexseq);

    
    print join("\t", "gene", "g_padj", "g_log2fc", "s_padj", "s_log2fc") . "\n";
    foreach my $gene (keys %genome_signif) {
        my $genome_best = $genome_signif{$gene};
        my $superT_best = $superT_signif{$gene};

        print join("\t", $gene, 
                   $genome_best->{padj}, $genome_best->{log2fold},
                   $superT_best->{padj}, $superT_best->{log2fold}) . "\n";
    }

    exit(0);
}
 

####
sub parse_signif_scores {
    my ($results_file) = @_;

    my %gene_to_scores;
    
    open(my $fh, $results_file) or die $!;
    my $header = <$fh>;
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $gene = $x[1];
        my $log2fold = $x[9];
        my $padj = $x[6];

        if ($padj eq "NA") { next; }
        
        push (@{$gene_to_scores{$gene}}, { log2fold => $log2fold,
                                           padj => $padj,
              } );

    }
    close $fh;

    my %gene_to_best;
    foreach my $gene (keys %gene_to_scores) {
        my @scores = @{$gene_to_scores{$gene}};

        @scores = sort {$a->{padj}<=>$b->{padj}
                        ||
                            abs($b->{log2fold}) <=> abs($a->{log2fold}) } @scores;

        my $best = shift @scores;
        $gene_to_best{$gene} = $best;
    }

    return(%gene_to_best);
}

