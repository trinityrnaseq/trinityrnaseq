#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use FindBin;
use Data::Dumper;

my $usage = <<__EOUSAGE__;

#################################################################################### 
#
# Required:
#
#  --matrix <string>       matrix.normalized.FPKM
#
# Optional:
#
#  -P <float>             p-value cutoff for FDR  (default: 0.001)
# 
#  -C <float>             min abs(log2(a/b)) fold change (default: 2  (meaning 2^(2) or 4-fold).
#
#  --output <float>       prefix for output file (default: "diffExpr.P\${Pvalue}_C\${C})
#
#
#
#
# Misc:
#
#  --samples <string>                     sample-to-replicate mappings (provided to run_DE_analysis.pl)
#
#  --max_DE_genes_per_comparison <int>    extract only up to the top number of DE features within each pairwise comparison.
#                                         This is useful when you have massive numbers of DE features but still want to make
#                                         useful heatmaps and other plots with more manageable numbers of data points.
#
#  --order_columns_by_samples_file        instead of clustering samples or replicates hierarchically based on gene expression patterns,
#                                         order columns according to order in the --samples file.
#
#  --max_genes_clust <int>                default: 10000  (if more than that, heatmaps are not generated, since too time consuming)
#
#  --examine_GO_enrichment                run GO enrichment analysis
#       --GO_annots <string>              GO annotations file
#       --gene_lengths <string>           lengths of genes file
#
#
##############################################################



__EOUSAGE__

    ;


my $matrix_file;
my $p_value = 0.001;
my $log2_fold_change = 2;
my $output_prefix = "";
my $FORCE_FLAG = 0;
my $help_flag = 0;
my $DESeq_mode = 0;

my $max_DE_genes_per_comparison;
my $max_genes_clust = 10000;

my $samples_file;

my $include_heatmaps_for_pairwise_comparisons = 0;

my $order_columns_by_samples_file = 0;

####
my $examine_GO_enrichment_flag;
my $GO_annots_file;
my $gene_lengths_file;



&GetOptions (  'h' => \$help_flag,
               
               'matrix=s' => \$matrix_file,
               'P=f' => \$p_value,
               'C=f' => \$log2_fold_change,
               'output=s' => \$output_prefix,
               'FORCE' => \$FORCE_FLAG, # for exploratory purposes.
               
               'max_DE_genes_per_comparison=i' => \$max_DE_genes_per_comparison,
               
               'max_genes_clust=i' => \$max_genes_clust,

               # gene ontology enrichment
               'examine_GO_enrichment' => \$examine_GO_enrichment_flag,
               'GO_annots=s' => \$GO_annots_file,
               'gene_lengths=s' => \$gene_lengths_file,
                 
               'samples=s' => \$samples_file,
             
               "order_columns_by_samples_file" => \$order_columns_by_samples_file,

               );


if ($help_flag) {
    die $usage;
}


unless ($matrix_file && -s $matrix_file) {
    die $usage;
}

if ($examine_GO_enrichment_flag && ! ($GO_annots_file && $gene_lengths_file)) {
    die "Error, need --GO_annots and --gene_lengths parameters specified for GO enrichment analysis";
}

unless ($output_prefix) {
    $output_prefix = "diffExpr.P${p_value}_C${log2_fold_change}";
}



main: {

    my @DE_result_files = <*.DE_results>;
    unless (@DE_result_files) {
        die "Error, no DE_results files!  This needs to be run in the edgeR or DESeq output directory";
    }

    my %column_header_to_index = &parse_column_headers($DE_result_files[0]);
    
    ## want P-value and log2(FC) columns
    
    my $pvalue_index = $column_header_to_index{padj}
    || $column_header_to_index{FDR}
    || die "Error, cannot identify FDR column from " . Dumper(\%column_header_to_index);
    

    my $log2FC_index = $column_header_to_index{log2FoldChange} 
    || $column_header_to_index{logFC} 
    || die "Error, cannot identify logFC column from " . Dumper(\%column_header_to_index);
    
    my $Mvalue_index = $column_header_to_index{baseMean}
    || $column_header_to_index{logCPM}
    || die "Error, cannot identify average counts column";
    
    
    ## store matrix
    my %read_count_rows;
    my $matrix_header;
    
    {
        # not counts, but FPKM values! 
        open (my $fh, $matrix_file) or die "Error, cannot read file $matrix_file";
        $matrix_header = <$fh>;
        chomp $matrix_header;
        $matrix_header =~ s/^\s+//;
                
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = shift @x;
            $read_count_rows{$acc} = join("\t", @x);
        }
       
    }
    my %samples_to_replicates;
    if ($samples_file) {
        open (my $fh, $samples_file) or die "Error, cannot open file $samples_file";
        while (<$fh>) {
            chomp;
            unless (/\w/) { next; }
            my ($sample, $replicate, @rest) = split(/\s+/, $_);
            if ($sample && $replicate) {
                
                
            push (@{$samples_to_replicates{$sample}}, $replicate);
            }
            else{
                print STDERR "-ignoring line of samples file: $_\n";
            }
        }
        close $fh;
    }
    else {
        foreach my $header_field (split(/\t/, $matrix_header)) {
            $samples_to_replicates{$header_field} = [$header_field];
        }
    }
    
    ## get list of genes that meet the criterion:
    my %diffExpr = &parse_result_files_find_diffExp(\@DE_result_files,
                                                    $p_value, $pvalue_index,
                                                    $log2_fold_change, $log2FC_index, \%read_count_rows, 
                                                    $matrix_header, $max_DE_genes_per_comparison, \%samples_to_replicates);
    
    unless (%diffExpr) {
        die "Error, no differentially expressed transcripts identified at cuttoffs: P:$p_value, C:$log2_fold_change";
    }

    my $diff_expr_matrix = "$output_prefix.matrix";
    if ($max_DE_genes_per_comparison) {
        $diff_expr_matrix = "$output_prefix.${max_DE_genes_per_comparison}_max_DE_eachcompare.matrix";
    }
    
    {
        open (my $ofh, ">$diff_expr_matrix") or die "Error, cannot write to file $diff_expr_matrix";
        print $ofh "\t$matrix_header\n";
        foreach my $acc (keys %diffExpr) {
            my $counts_row = $read_count_rows{$acc} or die "Error, no read counts row for $acc";
            print $ofh join("\t", $acc, $counts_row) . "\n";
        }
        close $ofh;
    }
    
    my $num_DE = scalar(keys %diffExpr);
    print STDERR "\n\n** Found $num_DE features as differentially expressed.\n\n";
    
    if ($num_DE > $max_genes_clust) {
        print STDERR "\n\n $num_DE exceeds $max_genes_clust so skipping clustering and heatmap construction.  Set --max_genes_clust to higher value to enable clustering, and excercise patience accordingly. :) \n\n\n";
    }
    else {
        &cluster_diff_expressed_transcripts($diff_expr_matrix);
    }


    exit(0);
    
}

####
sub cluster_diff_expressed_transcripts {
    my ($diff_expr_matrix_file) = @_;
    
    my $cmd = "$FindBin::Bin/PtR -m $diff_expr_matrix_file --log2 --heatmap --min_colSums 0 --min_rowSums 0 --gene_dist euclidean --sample_dist euclidean --sample_cor_matrix --center_rows --save @ARGV";
    
    if ($samples_file) {
        $cmd .= " -s $samples_file";
                
        if ($order_columns_by_samples_file) {
            $cmd .= " --order_columns_by_samples_file ";
        }

    }
    
    eval {
        &process_cmd($cmd);
    };

    if ($@) {
        print STDERR "Error, $@";
    }
    
    return;
    
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


####
sub parse_result_files_find_diffExp {
    my ($result_files_aref, 
        $max_p_value, $pvalue_index, 
        $min_abs_log2_fold_change, $log2FC_index, 
        $read_fpkm_rows_href, $fpkm_matrix_header, $max_DE_per_comparison, $samples_to_replicates_href) = @_;
    
    my %diff_expr;
    
    my %DE_pair_counts;
    
    foreach my $result_file (@$result_files_aref) {
        
        $result_file =~ /matrix\.(\S+)_vs_(\S+)\.[^\.]+\.DE_results$/ or die "Error, cannot extract condition names from $result_file";
        my ($condA, $condB) = ($1, $2);
        
        my $pairwise_samples_file = "$result_file.samples";
        { 
            ## write a samples file
            open (my $ofh, ">$pairwise_samples_file") or die "Error, cannot write to $pairwise_samples_file";
            foreach my $cond ($condA, $condB) {
                if (! exists $samples_to_replicates_href->{$cond}) {
                    die "\n\n\tError, not finding replicates corresponding to cond: $cond.  Please use the --samples parameter to specify sample/replicate relationships\n\n" . Dumper($samples_to_replicates_href);
                }
                my @replicates = @{$samples_to_replicates_href->{$cond}};
                foreach my $replicate (@replicates) {
                    print $ofh join("\t", $cond, $replicate) . "\n";
                }
            }
            close $ofh;
        }
        
        open (my $fh, $result_file) or die "Error, cannot open file $result_file";
        my $result_file_out_prefix = "$result_file.P${max_p_value}_C${min_abs_log2_fold_change}";
        
        my $condA_up_subset_file = "$result_file_out_prefix.$condA-UP.subset";
        my $condB_up_subset_file = "$result_file_out_prefix.$condB-UP.subset";

        open (my $condA_up_ofh, ">$condA_up_subset_file") or die "Error, cannot write to $condA_up_subset_file";
        open (my $condB_up_ofh, ">$result_file_out_prefix.$condB-UP.subset") or die "Error, cannot write to $condB_up_subset_file";
        
        my $count = 0;

        my $header = <$fh>;
        chomp $header;
        unless ($header =~ /^id\t/) {
            print $condA_up_ofh "id\t";
            print $condB_up_ofh "id\t";
        }
        print $condA_up_ofh "$header\t$fpkm_matrix_header\n";
        print $condB_up_ofh "$header\t$fpkm_matrix_header\n";

        my $countA = 0;
        my $countB = 0;

        my %condA_up_genes;
        my %condB_up_genes;
        
        my $rank = 0;
        while (<$fh>) {
            if (/^\#/) { next; }
            chomp;
            my $line = $_;
                        
            my @x = split(/\t/);
            my $log_fold_change = $x[$log2FC_index];
            my $fdr = $x[$pvalue_index];
            my $id = $x[0];
            
            if ($log_fold_change eq "NA") { next; }
            
            $rank++;
            
            if ( ($log_fold_change =~ /inf/i || abs($log_fold_change) >= $min_abs_log2_fold_change)
                 &&
                 $fdr <= $max_p_value) {
                
                $count++;
                if ((! $max_DE_per_comparison) || ($max_DE_per_comparison &&  $count <= $max_DE_per_comparison)) {
                    
                    $diff_expr{$id} = 1;
                    
                    my $matrix_counts = $read_fpkm_rows_href->{$id} || die "Error, no counts from matrix for $id";
            
                    if ($log_fold_change > 0) {
                        
                        print $condB_up_ofh "$line\t$matrix_counts\n";
                        $countB++;
                        $condB_up_genes{$id} = $rank;
                    }
                    else {
                        print $condA_up_ofh "$line\t$matrix_counts\n";
                        $countA++;
                        $condA_up_genes{$id} = $rank;
                    }
                }
            }
        }
        close $fh;
        close $condA_up_ofh;
        close $condB_up_ofh;
     
        ## generate heatmaps for the enriched genes
        if ($include_heatmaps_for_pairwise_comparisons) 
        {
            ## conditionA enriched heatmaps.
            &write_matrix_generate_heatmap("$condA_up_subset_file.matrix", \%condA_up_genes, $read_fpkm_rows_href, $fpkm_matrix_header, $pairwise_samples_file);
            
            # conditionB enriched heatmaps
            &write_matrix_generate_heatmap("$condB_up_subset_file.matrix", \%condB_up_genes, $read_fpkm_rows_href, $fpkm_matrix_header, $pairwise_samples_file);
            
        }
        
        ## do GO enrichment analysis
        if ($examine_GO_enrichment_flag) {
            
            my $cmd = "$FindBin::Bin/run_GOseq.pl --GO_assignments $GO_annots_file --lengths $gene_lengths_file --genes_single_factor $condA_up_subset_file";
            &process_cmd($cmd) if $countA;

            $cmd = "$FindBin::Bin/run_GOseq.pl --GO_assignments $GO_annots_file --lengths $gene_lengths_file --genes_single_factor $condB_up_subset_file";
            &process_cmd($cmd) if $countB;
            
        }
        
   
        $DE_pair_counts{$condA}->{$condB} = $count;
        $DE_pair_counts{$condB}->{$condA} = $count;
    }
    
    {
        # write DE count matrix:
        my $num_DE_counts_outfile = "numDE_feature_counts.P${max_p_value}_C${min_abs_log2_fold_change}.matrix";
        open (my $ofh, ">$num_DE_counts_outfile") or die "Error, cannot write to file $num_DE_counts_outfile";
        my @conditions = sort keys %DE_pair_counts;
        print $ofh "\t" . join("\t", @conditions) . "\n";
        foreach my $conditionA (@conditions) {
            print $ofh "$conditionA";
            foreach my $conditionB (@conditions) {
                my $count_DE = $DE_pair_counts{$conditionA}->{$conditionB} || 0;
                print $ofh "\t$count_DE";
            }
            print $ofh "\n";
        }
        close $ofh;
    }
    
    return(%diff_expr);
    
}


####
sub parse_column_headers {
    my ($DE_result_file) = @_;

    open (my $fh, $DE_result_file) or die "Error, cannot open file $DE_result_file";
    my $top_line = <$fh>;
    my $second_line = <$fh>;
    close $fh;

    chomp $top_line;
    my @columns = split(/\t/, $top_line);
    
    chomp $second_line;
    my @second_line_columns = split(/\t/, $second_line);
    if (scalar(@columns) == scalar(@second_line_columns) -1) {
        # weird R thing where the header can be off by one due to row.names
        unshift (@columns, "id");
    }
    
    my %indices;
    for (my $i = 0; $i <= $#columns; $i++) {
        $indices{$columns[$i]} = $i;
    }

    return(%indices);
}
    
####
sub write_matrix_generate_heatmap {
    my ($matrix_out_file, $genes_href, $fpkm_matrix_rows_href, $fpkm_matrix_header, $pairwise_samples_file) = @_;
    
    open (my $ofh, ">$matrix_out_file") or die "Error, cannot write to $matrix_out_file";
    print $ofh "$fpkm_matrix_header\n";
    
    foreach my $gene (sort {$genes_href->{$a}<=>$genes_href->{$b}} keys %$genes_href) { # output in rank of DE
        
        my $row = $fpkm_matrix_rows_href->{$gene} or die "Error, no fpkm row for gene $gene";
        print $ofh "$gene\t$row\n";
    }
    close $ofh;

    my $cmd = "$FindBin::Bin/PtR -m $matrix_out_file -s $pairwise_samples_file --log2 --heatmap --min_colSums 0 --min_rowSums 0 --gene_dist euclidean --sample_dist euclidean @ARGV ";
    
    if ($samples_file) {
        $cmd .= " -s $samples_file ";
        
        if ($order_columns_by_samples_file) {
            $cmd .= " --order_columns_by_samples_file ";
        }
        
    }
    
    

    &process_cmd($cmd);
 
    return;
}
