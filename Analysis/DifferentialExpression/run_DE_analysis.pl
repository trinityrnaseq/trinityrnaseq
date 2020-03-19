#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Cwd;
use FindBin;
use File::Basename;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;
use Data::Dumper;


my $ROTS_B = 500;
my $ROTS_K = 5000;

my $MIN_REPS_MIN_CPM = "2,1";


my $usage = <<__EOUSAGE__;


#################################################################################################
#
#  Required:
#
#  --matrix|m <string>               matrix of raw read counts (not normalized!)
#
#  --method <string>               edgeR|DESeq2|voom
#                                     note: you should have biological replicates.
#                                           edgeR will support having no bio replicates with
#                                           a fixed dispersion setting. 
#
#  Optional:
#
#  --samples_file|s <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1
#                                        cond_A    cond_A_rep2
#                                        cond_B    cond_B_rep1
#                                        cond_B    cond_B_rep2
#
#
#  General options:
#
#  --min_reps_min_cpm  <string>    default: $MIN_REPS_MIN_CPM  (format: 'min_reps,min_cpm')
#                                  At least min count of replicates must have cpm values > min cpm value.
#                                     (ie. filtMatrix = matrix[rowSums(cpm(matrix)> min_cpm) >= min_reps, ]  adapted from edgeR manual)
#                                      Note, ** if no --samples_file, default for min_reps is set = 1 **
#
#  --output|o                      name of directory to place outputs (default: \$method.\$pid.dir)
#
#  --reference_sample <string>     name of a sample to which all other samples should be compared.
#                                   (default is doing all pairwise-comparisons among samples)
#
#  --contrasts <string>            file (tab-delimited) containing the pairs of sample comparisons to perform.
#                                  ex. 
#                                       cond_A    cond_B
#                                       cond_Y    cond_Z
#
#
###############################################################################################
#
#  ## EdgeR-related parameters
#  ## (no biological replicates)
#
#  --dispersion <float>            edgeR dispersion value (Read edgeR manual to guide your value choice)
#                                    http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
#
###############################################################################################
#
#   Documentation and manuals for various DE methods.  Please read for more advanced and more
#   fine-tuned DE analysis than provided by this helper script.
#
#  edgeR:       http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
#  DESeq2:      http://bioconductor.org/packages/release/bioc/html/DESeq2.html    
#  voom/limma:  http://bioconductor.org/packages/release/bioc/html/limma.html
#
###############################################################################################



__EOUSAGE__


    ;


#  ## ROTS parameters
#  --ROTS_B <int>                   : number of bootstraps and permutation resampling (default: $ROTS_B)
#  --ROTS_K <int>                   : largest top genes size (default: $ROTS_K)
#
#  ROTS:        http://www.btk.fi/research/research-groups/elo/software/rots/


my $matrix_file;
my $method;
my $samples_file;

my $help_flag;
my $output_dir;
my $dispersion; # I've used 0.1 myself - but read the manual to guide your choice.

my $contrasts_file;

my $reference_sample;

my $make_tar_gz_file = 0;


&GetOptions ( 'h' => \$help_flag,
              'matrix|m=s' => \$matrix_file,              
              'method=s' => \$method,
              'samples_file|s=s' => \$samples_file,
              'output|o=s' => \$output_dir,
              'min_reps_min_cpm=s' => \$MIN_REPS_MIN_CPM,
              'dispersion=f' => \$dispersion,
              
              'reference_sample=s' => \$reference_sample,
              'contrasts=s' => \$contrasts_file,
              
              'tar_gz_outdir' => \$make_tar_gz_file,


              'ROTS_B=i' => \$ROTS_B,
              'ROTS_K=i' => \$ROTS_K,
    );



if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, don't understand options: @ARGV, please check spelling matches usage info.";
}


unless ($matrix_file 
        && $method
    ) { 
    
    die $usage;
    
}

if ($matrix_file =~ /fpkm|tpm/i) {
    die "Error, be sure you're using a matrix file that corresponds to raw counts, and not FPKM values.\n"
        . "If this is correct, then please rename your file, and remove fpkm or tpm from the name.\n\n";
}


unless ($method =~ /^(edgeR|DESeq2|voom|ROTS|GLM)$/) {
    die "Error, do not recognize --method [$method]";
}

my ($MIN_REPS, $MIN_CPM) = split(/,/, $MIN_REPS_MIN_CPM);

if ($samples_file) {
    unless ($MIN_REPS > 0 && $MIN_CPM > 0) {
        die "Error, --min_reps_min_cpm $MIN_REPS_MIN_CPM must include values > 0 in comma-delimited format. ex.  '2,1' ";
    }
}
else {
    print STDERR "-note, no biological replicates identified, so setting min reps = $MIN_REPS.\n";
    $MIN_REPS = 1;
}


main: {


    my $workdir = cwd();
    
    
    my %sample_name_to_column = &get_sample_name_to_column_index($matrix_file);
    
    my %samples;
    if ($samples_file) {
        unless ($samples_file =~ /^\//) {
            $samples_file = cwd() . "/$samples_file";
        }
        
        %samples = &parse_sample_info($samples_file);
    }
    else {
        # no replicates, so assign each sample to itself as a single replicate
        foreach my $sample_name (keys %sample_name_to_column) {
            $samples{$sample_name} = [$sample_name];
        }
    }

    print Dumper(\%samples);
        
    if ($matrix_file !~ /^\//) {
        ## make full path
        $matrix_file = cwd() . "/$matrix_file";
    }
        
    unless ($output_dir) {
        $output_dir = "$method.$$.dir";
    }
    
    unless (-d $output_dir) {
        mkdir($output_dir) or die "Error, cannot mkdir $output_dir";
    }
    chdir $output_dir or die "Error, cannot cd to $output_dir";
    

    my @sample_names = keys %samples;


    if ($method eq "GLM") {
        unless ($samples_file) { 
            die "Error, need samples file for GLM";
        }
        ## samples file here requires a different format:
        # replicate (tab) attrA [(tab) attrB, ...]

        &run_GLM($matrix_file, \%samples, \%sample_name_to_column);
    }
    else {
        # edgeR or DESeq pairwise comparison between samples:
        
        my @DE_contrasts;
        
        if ($reference_sample) {
            
            my @other_samples = grep { $_ ne $reference_sample} @sample_names;

            unless (@other_samples) {
                die "Error, couldn't extract non-reference samples from list: @sample_names";
            }

            foreach my $other_sample (@other_samples) {
                push (@DE_contrasts, [$reference_sample, $other_sample]);
            }
            
        }
        elsif ($contrasts_file) {
            
            unless ($contrasts_file =~ /^\//) {
                $contrasts_file = "$workdir/$contrasts_file";
            }
            
            open (my $fh, $contrasts_file) or die "Error, cannot open file $contrasts_file";
            while (<$fh>) {
                chomp;
                unless (/\w/) { next; }
                if (/^\#/) { next; }
                my ($sampleA, $sampleB) = split(/\s+/);
                unless ($sampleA && $sampleB) {
                    die "Error, didn't read a pair of tab-delimited samples from $contrasts_file, line: $_";
                }
                push (@DE_contrasts, [$sampleA, $sampleB]);
            }
            close $fh;
        }
        else {
            ## performing all pairwise comparisons:
            
            @sample_names = sort @sample_names;
            for (my $i = 0; $i < $#sample_names; $i++) {
                for (my $j = $i + 1; $j <= $#sample_names; $j++) {

                    push (@DE_contrasts, [$sample_names[$i], $sample_names[$j]]);
                }
            }
        }
        
        print STDERR "Contrasts to perform are: " . Dumper(\@DE_contrasts);
        
        foreach my $DE_contrast (@DE_contrasts) {
                            
            my ($sample_a, $sample_b) = @$DE_contrast;
            
            if ($method eq "edgeR") {
                &run_edgeR_sample_pair($matrix_file, \%samples, \%sample_name_to_column, $sample_a, $sample_b);
                
            }
            elsif ($method eq "DESeq2") {
                &run_DESeq2_sample_pair($matrix_file, \%samples, \%sample_name_to_column, $sample_a, $sample_b);
            }
            elsif ($method eq 'voom') {
                &run_limma_voom_sample_pair($matrix_file, \%samples, \%sample_name_to_column, $sample_a, $sample_b);
            }
            elsif ($method eq 'ROTS') {
                &run_ROTS_sample_pair($matrix_file, \%samples, \%sample_name_to_column, $sample_a, $sample_b);
            }
        }
    }

    if ($make_tar_gz_file) {
        chdir $workdir or die "Error, cannot cd to $workdir";
        my $cmd = "tar -zcvf $output_dir.tar.gz $output_dir";
        &process_cmd($cmd);
    }
            

    exit(0);
}

####
sub parse_sample_info {
    my ($sample_file) = @_;

    my %samples;

    open (my $fh, $sample_file) or die "Error, cannot locate samples file: $sample_file";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; } # allow comments
        chomp;
        s/^\s+//; # trim any leading ws
        my @x = split(/\s+/); # now ws instead of just tabs
        if (scalar @x < 2) { next; }
        my ($sample_name, $replicate_name, @rest) = @x;
        
        #$sample_name =~ s/^\s|\s+$//g;
        #$replicate_name =~ s/^\s|\s+$//g;
        
        push (@{$samples{$sample_name}}, $replicate_name);
    }
    close $fh;

    return(%samples);
}

####
sub get_sample_name_to_column_index {
    my ($matrix_file) = @_;

    my %column_index;

    open (my $fh, $matrix_file) or die "Error, cannot open file $matrix_file";
    my $header_line = <$fh>;

    $header_line =~ s/^\#//; # remove comment field.
    $header_line =~ s/^\s+|\s+$//g;
    my @samples = split(/\t/, $header_line);

    { # check for disconnect between header line and data lines
        my $next_line = <$fh>;
        my @x = split(/\t/, $next_line);
        print STDERR "Got " . scalar(@samples) . " samples, and got: " . scalar(@x) . " data fields.\n";
        print STDERR "Header: $header_line\nNext: $next_line\n";
        
        if (scalar(@x) == scalar(@samples)) {
            # problem... shift headers over, no need for gene column heading
            shift @samples;
            print STDERR "-shifting sample indices over.\n";
        }
    }
    close $fh;
            
    
    my $counter = 0;
    foreach my $sample (@samples) {
        $counter++;
        
        $sample =~ s/\.(isoforms|genes)\.results$//; 
        
        $column_index{$sample} = $counter;
    }

    use Data::Dumper;
    print STDERR Dumper(\%column_index);
    

    return(%column_index);
    
}


####
sub run_edgeR_sample_pair {
    my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B) = @_;
         
    my $output_prefix = basename($matrix_file) . "." . join("_vs_", ($sample_A, $sample_B));
        
    my $Rscript_name = "$output_prefix.$sample_A.vs.$sample_B.EdgeR.Rscript";
    
    my @reps_A = @{$samples_href->{$sample_A}};
    my @reps_B = @{$samples_href->{$sample_B}};

    my $num_rep_A = scalar(@reps_A);
    my $num_rep_B = scalar(@reps_B);
    
    my @rep_column_indices;
    foreach my $rep_name (@reps_A, @reps_B) {
        my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
        push (@rep_column_indices, $column_index);
    }
        

    ## write R-script to run edgeR
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";

    print $ofh "if (! require(edgeR)) {\n";
    print $ofh "   source(\"https://bioconductor.org/biocLite.R\")\n";
    print $ofh "   biocLite(\"edgeR\")\n";
    print $ofh "   library(edgeR)\n";
    print $ofh "}\n\n";
    
    
    print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
    print $ofh "col_ordering = c(" . join(",", @rep_column_indices) . ")\n";
    print $ofh "rnaseqMatrix = data[,col_ordering]\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > $MIN_CPM) >= $MIN_REPS,]\n";
    print $ofh "conditions = factor(c(rep(\"$sample_A\", $num_rep_A), rep(\"$sample_B\", $num_rep_B)))\n";
    print $ofh "\n";
    print $ofh "exp_study = DGEList(counts=rnaseqMatrix, group=conditions)\n";
    print $ofh "exp_study = calcNormFactors(exp_study)\n";
    
    if ($num_rep_A > 1 && $num_rep_B > 1) {
        #print $ofh "exp_study = estimateCommonDisp(exp_study)\n";
        #print $ofh "exp_study = estimateTagwiseDisp(exp_study)\n";
        print $ofh "exp_study = estimateDisp(exp_study)\n"; # new recommended way
        print $ofh "et = exactTest(exp_study, pair=c(\"$sample_A\", \"$sample_B\"))\n";
    }
    elsif (!$dispersion) {
	die "Error, cannot calculate dispersions due to lack of replicates. Specify a dispersion parameter --dispersion <float>. See help for details\n";
    }
    else {
        unless (defined $dispersion) {
            confess "Error, must set --dispersion <float> when using edgeR w/o bio replicates";

        }
        print $ofh "et = exactTest(exp_study, pair=c(\"$sample_A\", \"$sample_B\"), dispersion=$dispersion)\n";
    }
    print $ofh "tTags = topTags(et,n=NULL)\n";
    print $ofh "result_table = tTags\$table\n";
    print $ofh "result_table = data.frame(sampleA=\"$sample_A\", sampleB=\"$sample_B\", result_table)\n";
    
    ## reset logfc so it's A/B instead of B/A to be consistent with DESeq2
    print $ofh "result_table\$logFC = -1 * result_table\$logFC\n";
    
    print $ofh "write.table(result_table, file=\'$output_prefix.edgeR.DE_results\', sep='\t', quote=F, row.names=T)\n";
    print $ofh "write.table(rnaseqMatrix, file=\'$output_prefix.edgeR.count_matrix\', sep='\t', quote=F, row.names=T)\n";
    
    ## generate MA and Volcano plots
    print $ofh "source(\"$FindBin::RealBin/R/rnaseq_plot_funcs.R\")\n";
    print $ofh "pdf(\"$output_prefix.edgeR.DE_results.MA_n_Volcano.pdf\")\n";

    print $ofh "plot_MA_and_Volcano(rownames(result_table), result_table\$logCPM, result_table\$logFC, result_table\$FDR)\n";
    print $ofh "dev.off()\n";
    
    close $ofh;

    ## Run R-script
    #my $cmd = "R --no-save --no-restore --no-site-file --no-init-file -q < $Rscript_name";
    my $cmd = "Rscript $Rscript_name";
    
    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        print STDERR "$@\n\n";
        print STDERR "\n\nWARNING: This EdgeR comparison failed...\n\n";
        ## if this is due to data paucity, such as in small sample data sets, then ignore for now.
    }
    

    return;
}
        
sub run_DESeq2_sample_pair {
    my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B) = @_;
         
    my $output_prefix = basename($matrix_file) . "." . join("_vs_", ($sample_A, $sample_B));
        
    my $Rscript_name = "$output_prefix.DESeq2.Rscript";
    
    my @reps_A = @{$samples_href->{$sample_A}};
    my @reps_B = @{$samples_href->{$sample_B}};

    my $num_rep_A = scalar(@reps_A);
    my $num_rep_B = scalar(@reps_B);
    

    if ($num_rep_A < 2 || $num_rep_B < 2) {
        print STDERR "DESeq2 only supported here with biological replicates for each condition. Skipping: $sample_A vs. $sample_B *** \n\n";
        return;
    }
    
    my @rep_column_indices;
    foreach my $rep_name (@reps_A, @reps_B) {
        my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
        push (@rep_column_indices, $column_index);
    }
    

    ## write R-script to run DESeq
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    print $ofh "if (! require(edgeR)) {\n";
    print $ofh "   source(\"https://bioconductor.org/biocLite.R\")\n";
    print $ofh "   biocLite(\"edgeR\")\n";
    print $ofh "   library(edgeR)\n";
    print $ofh "}\n\n";
    print $ofh "if (! require(DESeq2)) {\n";
    print $ofh "   source(\"https://bioconductor.org/biocLite.R\")\n";
    print $ofh "   biocLite(\"DESeq2\")\n";
    print $ofh "   library(DESeq2)\n";
    print $ofh "}\n\n";

    print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
    print $ofh "col_ordering = c(" . join(",", @rep_column_indices) . ")\n";
    print $ofh "rnaseqMatrix = data[,col_ordering]\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > $MIN_CPM) >= $MIN_REPS,]\n";
    print $ofh "conditions = data.frame(conditions=factor(c(rep(\"$sample_A\", $num_rep_A), rep(\"$sample_B\", $num_rep_B))))\n";
    print $ofh "rownames(conditions) = colnames(rnaseqMatrix)\n";
    print $ofh "ddsFullCountTable <- DESeqDataSetFromMatrix(\n"
             . "    countData = rnaseqMatrix,\n"
             . "    colData = conditions,\n"
             . "    design = ~ conditions)\n";
    print $ofh "dds = DESeq(ddsFullCountTable)\n";

    print $ofh "contrast=c(\"conditions\",\"$sample_A\",\"$sample_B\")\n";
    print $ofh "res = results(dds, contrast)\n";
    

    # adj from: Carsten Kuenne, thx!
    ##recreates baseMeanA and baseMeanB columns that are not created by default DESeq2 anymore
    print $ofh "baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)\$conditions == \"$sample_A\"])\n";
    print $ofh "baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)\$conditions == \"$sample_B\"])\n";
    print $ofh "res = cbind(baseMeanA, baseMeanB, as.data.frame(res))\n";
 
    ##adds an “id” column headline for column 0
    print $ofh "res = cbind(sampleA=\"$sample_A\", sampleB=\"$sample_B\", as.data.frame(res))\n";

    print $ofh "res\$padj[is.na(res\$padj)]  <- 1\n"; # Carsten Kuenne

    print $ofh "res = as.data.frame(res[order(res\$pvalue),])\n"; # rank by pvalue
    
    ## output results
    print $ofh "write.table(res, file=\'$output_prefix.DESeq2.DE_results\', sep='\t', quote=FALSE)\n";
    print $ofh "write.table(rnaseqMatrix, file=\'$output_prefix.DESeq2.count_matrix\', sep='\t', quote=FALSE)\n";
    
    
    ## generate MA and Volcano plots
    print $ofh "source(\"$FindBin::RealBin/R/rnaseq_plot_funcs.R\")\n";
    print $ofh "pdf(\"$output_prefix.DESeq2.DE_results.MA_n_Volcano.pdf\")\n";
    print $ofh "plot_MA_and_Volcano(rownames(res), log2(res\$baseMean+1), res\$log2FoldChange, res\$padj)\n";
    print $ofh "dev.off()\n";
        
    
    close $ofh;
    
    ## Run R-script
    #my $cmd = "R --no-save --no-restore --no-site-file --no-init-file -q < $Rscript_name";
    my $cmd = "Rscript $Rscript_name";
    
    &process_cmd($cmd);
    
    return;
}


####
sub run_limma_voom_sample_pair {
    my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B) = @_;
    
    my $output_prefix = basename($matrix_file) . "." . join("_vs_", ($sample_A, $sample_B));
        
    my $Rscript_name = "$output_prefix.$sample_A.vs.$sample_B.voom.Rscript";
    
    my @reps_A = @{$samples_href->{$sample_A}};
    my @reps_B = @{$samples_href->{$sample_B}};

    my $num_rep_A = scalar(@reps_A);
    my $num_rep_B = scalar(@reps_B);
    
    unless ($num_rep_A > 1 && $num_rep_B > 1) {
        die "Error, need multiple biological replicates for each sample in order to run voom";
    }

    my @rep_column_indices;
    foreach my $rep_name (@reps_A, @reps_B) {
        my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
        push (@rep_column_indices, $column_index);
    }
        

    ## write R-script to run edgeR
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    
    print $ofh "library(edgeR)\n";
    print $ofh "library(limma)\n";
    
    print $ofh "\n";
    
    print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
    print $ofh "col_ordering = c(" . join(",", @rep_column_indices) . ")\n";
    print $ofh "rnaseqMatrix = data[,col_ordering]\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > $MIN_CPM) >= $MIN_REPS,]\n";
    print $ofh "conditions = factor(c(rep(\"$sample_A\", $num_rep_A), rep(\"$sample_B\", $num_rep_B)))\n";
    print $ofh "\n";
    print $ofh "design = model.matrix(~conditions)\n";
    print $ofh "## TMM normalize data\n";
    print $ofh "lib_sizes = colSums(rnaseqMatrix)\n";
    print $ofh "tmm_norm_factors = calcNormFactors(rnaseqMatrix, method='TMM')\n";
    print $ofh "x = DGEList(counts=rnaseqMatrix)\n";
    print $ofh "# voom transformation\n";
    print $ofh "y = voom(x, design, lib.size=lib_sizes*tmm_norm_factors, plot=F)\n";
    print $ofh "fit = eBayes(lmFit(y,design))\n";
    print $ofh "tTags = topTable(fit,coef=2,number=Inf)\n";
    print $ofh "# output results, including average expression val for each feature\n";
    print $ofh "c = cpm(x)\n";
    print $ofh "m = apply(c, 1, mean)\n";
    print $ofh "tTags\$logFC = -1 * tTags\$logFC  # make A/B instead of B/A\n";
    print $ofh "tTags2 = cbind(tTags, logCPM=log2(m[rownames(tTags)]))\n";
    print $ofh "DE_matrix = data.frame(sampleA=\"$sample_A\", sampleB=\"$sample_B\", logFC=tTags\$logFC, logCPM=tTags2\$logCPM, PValue=tTags\$'P.Value', FDR=tTags\$'adj.P.Val')\n";
    print $ofh "rownames(DE_matrix) = rownames(tTags)\n";
    print $ofh "write.table(DE_matrix, file=\'$output_prefix.voom.DE_results\', sep='\t', quote=F, row.names=T)\n";
    print $ofh "write.table(rnaseqMatrix, file=\'$output_prefix.voom.count_matrix\', sep='\t', quote=F, row.names=T)\n";
    
    ## generate MA and Volcano plots
    print $ofh "# MA and volcano plots\n";
    print $ofh "source(\"$FindBin::RealBin/R/rnaseq_plot_funcs.R\")\n";
    print $ofh "pdf(\"$output_prefix.voom.DE_results.MA_n_Volcano.pdf\")\n";
    print $ofh "plot_MA_and_Volcano(rownames(tTags2), tTags2\$logCPM, tTags\$logFC, tTags\$'adj.P.Val')\n";
    print $ofh "dev.off()\n";
    
    close $ofh;

    ## Run R-script
    #my $cmd = "R --no-save --no-restore --no-site-file --no-init-file -q < $Rscript_name";
    my $cmd = "Rscript $Rscript_name";
    
    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        print STDERR "$@\n\n";
        print STDERR "\n\nWARNING: This voom comparison failed...\n\n";
        ## if this is due to data paucity, such as in small sample data sets, then ignore for now.
    }
    

    return;
}


####
sub run_ROTS_sample_pair {
    my ($matrix_file, $samples_href, $sample_name_to_column_index_href, $sample_A, $sample_B) = @_;
    
    my $output_prefix = basename($matrix_file) . "." . join("_vs_", ($sample_A, $sample_B));
        
    my $Rscript_name = "$output_prefix.$sample_A.vs.$sample_B.ROTS.Rscript";
    
    my @reps_A = @{$samples_href->{$sample_A}};
    my @reps_B = @{$samples_href->{$sample_B}};

    my $num_rep_A = scalar(@reps_A);
    my $num_rep_B = scalar(@reps_B);
    
    unless ($num_rep_A > 1 && $num_rep_B > 1) {
        die "Error, need multiple biological replicates for each sample in order to run ROTS";
    }

    my @rep_column_indices;
    foreach my $rep_name (@reps_A, @reps_B) {
        my $column_index = $sample_name_to_column_index_href->{$rep_name} or die "Error, cannot determine column index for replicate name [$rep_name]" . Dumper($sample_name_to_column_index_href);
        push (@rep_column_indices, $column_index);
    }
        

    ## write R-script to run DESeq2
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    
    print $ofh "library(edgeR)\n";
    print $ofh "library(limma)\n";
    print $ofh "library(ROTS)\n";
    
    print $ofh "\n";
    
    print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
    print $ofh "col_ordering = c(" . join(",", @rep_column_indices) . ")\n";
    print $ofh "rnaseqMatrix = data[,col_ordering]\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > $MIN_CPM) >= $MIN_REPS,]\n";
    print $ofh "conditions = factor(c(rep(\"$sample_A\", $num_rep_A), rep(\"$sample_B\", $num_rep_B)))\n";
    print $ofh "\n";
    print $ofh "design = model.matrix(~conditions)\n";
    print $ofh "## TMM normalize data\n";
    print $ofh "lib_sizes = colSums(rnaseqMatrix)\n";
    print $ofh "tmm_norm_factors = calcNormFactors(rnaseqMatrix, method='TMM')\n";
    print $ofh "x = DGEList(counts=rnaseqMatrix)\n";
    print $ofh "# voom transformation and ROTS (code derived from ROTS paper supp R code)\n";
    print $ofh "voom.data = voom(x, design, lib.size=lib_sizes*tmm_norm_factors, plot=F)\n";
    print $ofh "input_voom = voom.data\$E\n";
    print $ofh "# run ROTS for DE analysis\n";
    print $ofh "res_voom <- ROTS(data=input_voom,groups=as.numeric(conditions),B=$ROTS_B, K=$ROTS_K)\n";
    print $ofh "results = summary(res_voom, fdr=0.1)\n";

    print $ofh "# add logFC and logCPM to result table.\n";
    print $ofh "c = cpm(x)\n";
    print $ofh "m = apply(c, 1, mean)\n";
    print $ofh "sampleA_cpm_matrix = c[,conditions \%in% \"$sample_A\"]\n";
    print $ofh "mean_sampleA_cpm = apply(sampleA_cpm_matrix, 1, mean)\n";
    print $ofh "sampleB_cpm_matrix = c[,conditions \%in% \"$sample_B\"]\n";
    print $ofh "mean_sampleB_cpm = apply(sampleB_cpm_matrix, 1, mean)\n";
    print $ofh "pseudocount_cpm = 1\n";
    print $ofh "FC = (mean_sampleA_cpm + pseudocount_cpm) / (mean_sampleB_cpm + pseudocount_cpm)\n";
    print $ofh "logFC = log2(FC)\n";
    print $ofh "results = summary(res_voom, fdr=0.1)\n";
    print $ofh "feature_order = rownames(results)\n";
    print $ofh "final_table = data.frame(sampleA=\"$sample_A\", sampleB=\"$sample_B\", logCPM=log2(m+1)[feature_order], CPM_A=mean_sampleA_cpm[feature_order], CPM_B=mean_sampleB_cpm[feature_order], logFC=logFC[feature_order], results)\n";
    
    print $ofh "write.table(final_table, file=\"$output_prefix.ROTS.DE_results\", quote=F, sep='\t')\n";
    print $ofh "write.table(rnaseqMatrix, file=\"$output_prefix.ROTS.count_matrix\", quote=F, sep='\t')\n";
    
    
    ## generate MA and Volcano plots
    print $ofh "# MA and volcano plots\n";
    print $ofh "source(\"$FindBin::RealBin/R/rnaseq_plot_funcs.R\")\n";
    print $ofh "pdf(\"$output_prefix.voom.DE_results.MA_n_Volcano.pdf\")\n";
    print $ofh "plot_MA_and_Volcano(rownames(final_table), final_table\$logCPM, final_table\$logFC, final_table\$FDR)\n";
    print $ofh "dev.off()\n";
    
    
    close $ofh;

    ## Run R-script
    #my $cmd = "R --no-save --no-restore --no-site-file --no-init-file -q < $Rscript_name";
    my $cmd = "Rscript $Rscript_name";
    
    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        print STDERR "$@\n\n";
        print STDERR "\n\nWARNING: This ROTS comparison failed...\n\n";
        ## if this is due to data paucity, such as in small sample data sets, then ignore for now.
    }
    

    return;
}




####
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret ($ret) ";
    }

    return;
}

####
sub run_GLM {
    my ($matrix_file, $samples_href, $sample_name_to_column_index_href) = @_;
    

    my $output_prefix = basename($matrix_file);
                 
    my $Rscript_name = "$output_prefix.GLM.Rscript";
    
    ## write R-script to run edgeR
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    
    print $ofh "library(edgeR)\n";
    
    print $ofh "\n";
    
    print $ofh "design_matrix = read.table(\"$samples_file\", header=T, row.names=1)\n";
    print $ofh "groups = factor(apply(design_matrix, 1, paste, collapse='.'))\n";
    print $ofh "design_matrix = cbind(design_matrix, groups=groups)\n";
    
    print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
    print $ofh "rnaseqMatrix = round(data)\n";
    print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > $MIN_CPM) >= $MIN_REPS,]\n";


    print $ofh "design = model.matrix(~0+groups)\n";
    print $ofh "colnames(design) = levels(groups)\n";
    print $ofh "rownames(design) = rownames(design_matrix)\n";
    print $ofh "rnaseqMatrix = rnaseqMatrix[,rownames(design)] # ensure properly ordered according to design\n";
    print $ofh "exp_study = DGEList(counts=rnaseqMatrix, group=groups)\n";
    print $ofh "exp_study = estimateGLMCommonDisp(exp_study,design)\n";
    print $ofh "exp_study = estimateGLMTrendedDisp(exp_study, design)\n";
    print $ofh "exp_study = estimateGLMTagwiseDisp(exp_study, design)\n";
    print $ofh "fit = glmFit(exp_study, design)\n";
    print $ofh "## define your contrasts:\n";
    print $ofh "levels(groups) # examine the factor combinations\n";
    print $ofh "# contrast = makeContrasts((wt.T15-wt.T0)-(zcf15.T15-zcf15.T0), levels=design)\n";
    print $ofh "# lrt = glmLRT(fit, contrast=contrast)\n";
    print $ofh "# topTags(lrt, n=100)\n";
    
    
    close $ofh;

    
    return;
}
