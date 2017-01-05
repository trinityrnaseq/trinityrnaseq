#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

###############################################################################################
#
#  --factor_labeling <string>       tab delimited file with format:  factor<tab>feature_id
#   or
#  --genes_single_factor <string>   list of genes to test (can be a matrix, only the first column is used for gene IDs)
#
#  --GO_assignments <string>        extracted GO assignments with format: feature_id <tab> GO:000001,GO:00002,...
#
#  --lengths <string>               feature lengths file with format:  feature_id <tab> length
#
#  --background <string>            gene ids file that defines the full population of genes to consider in testing.
#                                   Ideally, these represent the genes that are expressed and relevant to this test
#                                   as opposed to using all genes in the genome.
#
###############################################################################################



__EOUSAGE__

    ;


my ($factor_labeling, $GO_file, $help_flag, $lengths_file, $genes_single_factor_file);
my $background_file;

&GetOptions("factor_labeling=s" => \$factor_labeling,
            "GO_assignments=s" => \$GO_file,
            "lengths=s" => \$lengths_file,
            
            "genes_single_factor=s" => \$genes_single_factor_file,
            "background=s" => \$background_file,
            "help|h" => \$help_flag,

    );


if ($help_flag) {
    die $usage;
}


unless (($factor_labeling || $genes_single_factor_file) && $GO_file && $lengths_file && $background_file) {
    die $usage;
}


main: {

    my $Rscript = "__runGOseq.R";
    open (my $ofh, ">$Rscript") or die $!;
    
    print $ofh "library(goseq)\n";
    print $ofh "library(GO.db)\n";
    print $ofh "library(qvalue)\n";
    


    print $ofh "# capture list of genes for functional enrichment testing\n";
    if ($genes_single_factor_file) {
        print $ofh "factor_labeling = read.table(\"$genes_single_factor_file\", row.names=1)\n";
        print $ofh "factor_labeling[,1] = rep('custom_list', dim(factor_labeling)[1])\n";
        print $ofh "factor_labeling = factor_labeling[,1,drop=F]\n";
    }
    else {
        print $ofh "factor_labeling = read.table(\"$factor_labeling\", row.names=2, header=F)\n";
    }
    
    print $ofh "colnames(factor_labeling) = c('type')\n";
    print $ofh "factor_list = unique(factor_labeling[,1])\n";
    print $ofh "DE_genes = rownames(factor_labeling)\n";
    

    print $ofh "\n\n# get gene lengths\n";
    print $ofh "gene_lengths = read.table(\"$lengths_file\", header=T, row.names=1, com='')\n";
    print $ofh "gene_lengths = as.matrix(gene_lengths[,1,drop=F])\n";

    print $ofh "\n\n# get background gene list\n";
    print $ofh "background = read.table(\"$background_file\", header=T, row.names=1)\n";
    print $ofh "background.gene_ids = rownames(background)\n";
    print $ofh "background.gene_ids = unique(c(background.gene_ids, DE_genes))\n"; 

    print $ofh "\n\n# parse GO assignments\n";
    print $ofh "GO_info = read.table(\"$GO_file\", header=F, row.names=1,stringsAsFactors=F)\n";
    
    print $ofh "GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))\n";
    print $ofh "names(GO_info_listed) = rownames(GO_info)\n";

    print $ofh "get_GO_term_descr =  function(x) {\n";
    print $ofh "    d = 'none';\n"
             . "    go_info = GOTERM[[x]];\n"
             . "    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}\n"
             . "    return(d);\n"
             . "}\n";


    
    print $ofh "\n\n#organize go_id -> list of genes\n";  
    print $ofh "GO_to_gene_list = list()\n";
    print $ofh "for (gene_id in names(GO_info_listed)) {\n";
    print $ofh "    go_list = GO_info_listed[[gene_id]]\n";
    print $ofh "    for (go_id in go_list) {\n";
    print $ofh "        GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)\n";
    print $ofh "    }\n";
    print $ofh "}\n";
    
    
    print $ofh "\n\n# GO-Seq protocol: build pwf based on ALL DE features\n";
    
    
    print $ofh "sample_set_gene_ids = background.gene_ids\n";
    print $ofh "sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]\n";
    print $ofh "GO_info_listed = GO_info_listed[ names(GO_info_listed) %in% sample_set_gene_ids ]\n";
    
    print $ofh "cat_genes_vec = as.integer(sample_set_gene_ids %in% rownames(factor_labeling))\n";
    
    print $ofh "pwf=nullp(cat_genes_vec, bias.data=sample_set_gene_lengths)\n";
    print $ofh "rownames(pwf) = sample_set_gene_ids\n";

    
    print $ofh "\n\n# perform functional enrichment testing for each category.\n";
    print $ofh "for (feature_cat in factor_list) {\n";
    print $ofh "   message('Processing category: ', feature_cat)\n";
    print $ofh "    gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling\$type == feature_cat]\n";
    print $ofh "    cat_genes_vec = as.integer(sample_set_gene_ids %in% gene_ids_in_feature_cat)\n";
    print $ofh "    pwf\$DEgenes = cat_genes_vec\n";
    print $ofh "    res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=TRUE)\n";
    
    ## Process the over-represented    
    print $ofh "    ## over-represented categories:\n";
    print $ofh "     pvals = res\$over_represented_pvalue\n";
    print $ofh "     pvals[pvals > 1 - 1e-10] = 1 - 1e-10\n";
    print $ofh "     q = qvalue(pvals)\n";
    print $ofh "     res\$over_represented_FDR = q\$qvalues\n";
    
    if ($genes_single_factor_file) {
        print $ofh "go_enrich_filename = paste(\"$genes_single_factor_file\", '.GOseq.enriched', sep='')\n";
    }
    else {
        print $ofh "    go_enrich_filename = paste(feature_cat,'.GOseq.enriched', sep='')\n";
    }
    print $ofh "    result_table = res[res\$over_represented_pvalue<=0.05,]\n";

    print $ofh "    descr = unlist(lapply(result_table\$category, get_GO_term_descr))\n";
    print $ofh "    result_table\$go_term = descr;\n";

    print $ofh "    result_table\$gene_ids = do.call(rbind, lapply(result_table\$category, function(x) { \n" .
               "            gene_list = GO_to_gene_list[[x]]\n" .
               "            gene_list = gene_list[gene_list %in% rownames(factor_labeling)]\n" .
               "            paste(gene_list, collapse=', ');\n" .
               "     }) )\n";
    

    print $ofh "    write.table(result_table[order(result_table\$over_represented_pvalue),], file=go_enrich_filename, sep='\t', quote=F, row.names=F)\n";
    

    ## Process the under-represented    
    print $ofh "    ## under-represented categories:\n";
    
    print $ofh "     pvals = res\$under_represented_pvalue\n";
    print $ofh "     pvals[pvals>1-1e-10] = 1 - 1e-10\n";
    print $ofh "     q = qvalue(pvals)\n";
    print $ofh "     res\$under_represented_FDR = q\$qvalues\n";
    
    if ($genes_single_factor_file) {
        print $ofh "    go_depleted_filename = paste(\"$genes_single_factor_file\", '.GOseq.depleted', sep='')\n";
    }
    else {
        print $ofh "    go_depleted_filename = paste(feature_cat,'.GOseq.depleted', sep='')\n";
    }
    
    print $ofh "    result_table = res[res\$under_represented_pvalue<=0.05,]\n";
    
    print $ofh "    descr = unlist(lapply(result_table\$category, get_GO_term_descr))\n";
    print $ofh "    result_table\$go_term = descr;\n";
    print $ofh "    write.table(result_table[order(result_table\$under_represented_pvalue),], file=go_depleted_filename, sep='\t', quote=F, row.names=F)\n";
    
            
    print $ofh "}\n";
    
    close $ofh;

    my $cmd = "R --no-save --no-restore --no-site-file --no-init-file --quiet < $Rscript";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    else {
        print STDERR "\n\nDone.\n\n";
    }
    

    exit(0);
}

__END__

  

Notes:

1. Get the transcript GO annotation by running Trinotate, getting a trinotate.xls report file, and then running:

     trinotate-code/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate.xls -G --include_ancestral_terms > go_annotations.txt

    # use -T instead of -G in above to get transcript instead of gene-level annotations.


2.  Run GO-Seq like so, using this script 'run_GOseq.pl' included in Trinity:

      TRINITY_HOME/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling  factor_labeling.txt  --GO_assignments go_annotations.txt --lengths gene.lengths.txt

       The 'factor_labeling.txt' file should be of format:

              gene_id (tab) factor

       where factor is a string describing that subset of genes.

       For example:

              my_gene_A (tab) diff_expressed_cond_X_Y
              my_gene_B (tab) diff_expressed_cond_X_Y
              ...
              my_gene_M (tab) diff_cond_W_Z
              my_gene_N (tab) diff_cond_W_Z
              ...

      You can come up with whatever gene subset you want and call it whatever you want.  The enrichment tests will be performed separately for 
      each factor defined.

      The gene.lengths.txt file has the format

       gene (tab) length

        and you can use the same file you used earlier as part of doing the TMM normalization step and generating your FPKM matrix.


               
