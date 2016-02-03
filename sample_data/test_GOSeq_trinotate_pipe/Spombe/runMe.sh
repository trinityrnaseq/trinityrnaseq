#!/bin/sh
../../../Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling hs_induced_vs_log.factors --GO_assignments Trinotate_report.xls.trans.gene_ontology --lengths Trinity.seq_lengths

../../../Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling ds_induced_vs_log.factors --GO_assignments Trinotate_report.xls.trans.gene_ontology --lengths Trinity.seq_lengths
