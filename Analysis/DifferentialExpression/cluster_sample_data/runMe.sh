#!/bin/bash -ve

PROG=`dirname $0`/../PtR

if [ -e MLF_ESC_NPC.cuff.genes.fpkm.matrix.gz ] && [ ! -e MLF_ESC_NPC.cuff.genes.fpkm.matrix ]; then
    gunzip -c MLF_ESC_NPC.cuff.genes.fpkm.matrix.gz > MLF_ESC_NPC.cuff.genes.fpkm.matrix
fi

# assessing read counts and genes mapped
$PROG --matrix MLF_ESC_NPC.cuff.genes.fpkm.matrix --samples samples.txt --boxplot_log2_dist 1 --output boxplots

# compare replicates
$PROG --matrix MLF_ESC_NPC.cuff.genes.fpkm.matrix --samples samples.txt --compare_replicates --log2  --output rep_compare

# correlation of expression across all samples
$PROG --matrix MLF_ESC_NPC.cuff.genes.fpkm.matrix --samples samples.txt --log2 --sample_cor_matrix --output sample_cor

# PCA analysis
$PROG --matrix MLF_ESC_NPC.cuff.genes.fpkm.matrix --samples samples.txt --log2 --prin_comp 4 --output prin_comp

# Most variably expressed genes
$PROG --matrix MLF_ESC_NPC.cuff.genes.fpkm.matrix --samples samples.txt --log2 --top_variable_genes 1000 --var_gene_method anova --heatmap --center_rows --output anovaTop1k

# combine analyses
$PROG --matrix MLF_ESC_NPC.cuff.genes.fpkm.matrix --samples samples.txt --log2 \
      --top_variable_genes 1000 --var_gene_method anova --output anovaTop1k --heatmap --prin_comp 3  --add_prin_comp_heatmaps 30 \
      --center_rows --output top1kvarPC3Gene30


