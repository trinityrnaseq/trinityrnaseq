
library("ctc");
library("ape");

get_cluster_info = function( Rdata_file ) {

    load(Rdata_file);

    ordered_genes_file = paste(Rdata_file, ".ordered_gene_matrix", sep='');
    ordered_genes = hc_genes$data;
    write.table(ordered_genes, file=ordered_genes_file, quote=F, sep="\t");
    
    gene_tree = hc2Newick(hc_genes);
    gene_tree_filename = paste(Rdata_file, ".gene_tree", sep='');
    write(gene_tree, file=gene_tree_filename);

    # get rid of the distances since these can sometimes cause problems with other software tools. 
    gene_nodist_tree_filename = paste(Rdata_file, ".gene_nodist_tree", sep='');
    t = read.tree(text=gene_tree);
    t$edge.length = NULL;
    write.tree(t, file=gene_nodist_tree_filename);
        
    sample_tree = hc2Newick(hc_samples);
    sample_tree_filename = paste(Rdata_file, ".sample_tree", sep='');
    write(sample_tree, file=sample_tree_filename);

    
}    

       

