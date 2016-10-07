#!/bin/bash

set -e


$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl -m Trinity_genes.counts.matrix -s samples.txt --method edgeR --dispersion 0.001 --contrasts contrasts.txt -o edgeR_out && \
    cd edgeR_out && \
    $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl  -m ../Trinity_genes.TMM.EXPR.matrix -s ../samples.txt --examine_GO_enrichment --GO_annots ../Trinotate_report.xls.gene_ontology --gene_lengths ../Trinity.gene.lengths

