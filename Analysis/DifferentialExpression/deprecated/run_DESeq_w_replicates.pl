#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use File::Basename;

use Data::Dumper;


my $usage = <<_EOUSAGE_;

#############################################################################################################
#
# --counts_matrix <string>    fragment counts per feature (as estimated by RSEM)
#
# --repA_name <string>        a name for condition A (anything will do)
# --repA_list <string>        comma-delmited list of column names corresponding to replicates of condition A
#
# --repB_name <string>        a name for condition B
# --repB_list <string>        comma-delmited list of column names corresponding to replicates of condition B
#
#
#
##############################################################################################################



_EOUSAGE_

    ;

my $counts_matrix_file;

my $repA_name;
my $repA_list;

my $repB_name;
my $repB_list;


&GetOptions ( 'counts_matrix=s' => \$counts_matrix_file,
              
              'repA_name=s' => \$repA_name,
              'repA_list=s' => \$repA_list,

              'repB_name=s' => \$repB_name,
              'repB_list=s' => \$repB_list,
              
              );


unless ($counts_matrix_file 
        && $repA_name && $repB_name
        && $repA_list && $repB_list) {
    die $usage;
}

unless ($repA_list =~ /,/) {
    die "Error, repA_list requires multiple column names";
}
unless ($repB_list =~ /,/) {
    die "Error, repB_list requires multiple column names";
}


main: {

    my @repA = split(/,/, $repA_list);
    my @repB = split(/,/, $repB_list);
    foreach my $rep (@repA, @repB) {
        $rep =~ s/\s+//g;
    }
    
    ## figure out column headers
    open (my $fh, $counts_matrix_file) or die "Error, cannot open file $counts_matrix_file";
    my $header = <$fh>;
    chomp $header;
    my @x = split(/\t/, $header);
    my %col_name_to_index;
    for (my $i = 1; $i <= $#x; $i++) {
        my $col_name = $x[$i];
        if (exists ($col_name_to_index{$col_name})) {
            die "Error, column headings must be unique";
        }

        $col_name_to_index{$col_name} = $i;
    }
    close $fh;
    
    
    ## verify column headings
    my @rep_indices;
    foreach my $rep (@repA, @repB) {
        if (my $col_no = $col_name_to_index{$rep}) {
            push (@rep_indices, $col_no);
        }
        else {
            die "Error, cannot locate heading \"$rep\" as a column header.  Headers are: " . join(" ", keys %col_name_to_index) . " ";
        }
    }
    
    
    my $output_prefix = basename($counts_matrix_file) . "." . join("_vs_", sort ($repA_name, $repB_name));
    
    my $Rscript_name = "$output_prefix.DESeq.Rscript";
    

    my $num_rep_A = scalar(@repA);
    my $num_rep_B = scalar(@repB);

    ## write R-script to run DESeq
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    
    print $ofh "library(DESeq)\n";
    print $ofh "library(qvalue)\n";
    print $ofh "\n";
    
    ## write MA-plot function (got all this code from Zehua)
    print $ofh "# function: MA plot to show differential genes by a given q value cutoff\n";
    print $ofh "# genes meet the q value cutoff are in red;\n";
    print $ofh "# genes don't meet the q value cutoff but have log2(fold) greater than 2 are in blue;\n";
    print $ofh "# all other genes are in grey\n";
    
    # it's built-in now.
    print $ofh "plotMA <- function(res, qvalue) {\n";
    print $ofh "plot((log2(res\$baseMeanA)+log2(res\$baseMeanB))/2, res\$log2FoldChange, pch=20, cex=.4, xlab=\"A\", ylab=\"M\", col=ifelse(res\$padj<qvalue & (res\$log2FoldChange > 1 | res\$log2FoldChange < -1), \"red\", ifelse(res\$log2FoldChange > 2 | res\$log2FoldChange < -2, \"blue\", \"dark grey\")), main=title, cex.main=1.5, cex.lab=1.2)\n";
    
    print $ofh "points(c(-5,20), c(1,1), type=\"l\", col=\"dark grey\")\n";
    print $ofh "points(c(-5,20), c(-1,-1), type=\"l\", col=\"dark grey\")\n";
    print $ofh "points(c(-5,20), c(2,2), type=\"l\", col=\"dark grey\")\n";
    print $ofh "points(c(-5,20), c(-2,-2), type=\"l\", col=\"dark grey\")\n";
    print $ofh "}\n";
    
    print $ofh "data = read.table(\"$counts_matrix_file\", header=T, row.names=1)\n";
    print $ofh "col_ordering = c(" . join(",", @rep_indices) . ")\n";
    print $ofh "rnaseqMatrix = data[,col_ordering]\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "conditions = factor(c(rep(\"$repA_name\", $num_rep_A), rep(\"$repB_name\", $num_rep_B)))\n";
    print $ofh "\n";
    print $ofh "exp_study = newCountDataSet(rnaseqMatrix, conditions)\n";
    print $ofh "exp_study = estimateSizeFactors(exp_study)\n";
    print $ofh "sizeFactors(exp_study)\n";
    #print $ofh "exp_study = estimateVarianceFunctions(exp_study)\n";
    print $ofh "exp_study = estimateDispersions(exp_study, sharingMode=\"gene-est-only\")\n";
    print $ofh "str(fitInfo(exp_study))\n";
    #print $ofh "plotDispEsts(exp_study)\n";
    print $ofh "\n";
    print $ofh "res = nbinomTest(exp_study, \"$repA_name\", \"$repB_name\")\n";
    print $ofh "\n";
   
    print $ofh "# Adding Qvalues\n";
    print $ofh "p = res\$pval\n";
    print $ofh "p[is.na(p)]=1\n";
    print $ofh "q = qvalue(p=p)\n";
    print $ofh "res = cbind(res, q\$qvalues)\n";
    

    ## output results
    print $ofh "write.table(res[order(res\$pval),], file=\'$output_prefix.DESeq.output\', sep='\t', quote=FALSE, row.names=FALSE)\n";
    
    ## make MA plot
    print $ofh "postscript(file=\"$output_prefix.DESeq.MAplot.eps\", horizontal=FALSE, width=7, height=7, paper=\"special\")\n";
    print $ofh "title=\"$repA_name vs. $repB_name\"\n";
    print $ofh "qvalue=0.05\n";
    print $ofh "plotMA(res, qvalue)\n";
    print $ofh "dev.off()\n";
    
    close $ofh;

    ## Run R-script
    my $cmd = "R --vanilla -q < $Rscript_name";
    &process_cmd($cmd);
    
    exit(0);
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
    
    
