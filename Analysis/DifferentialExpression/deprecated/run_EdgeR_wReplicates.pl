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
    
    my $Rscript_name = "$output_prefix.EdgeR.Rscript";
    

    my $num_rep_A = scalar(@repA);
    my $num_rep_B = scalar(@repB);

    ## write R-script to run DESeq
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    
    print $ofh "library(edgeR)\n";

    print $ofh "\n";
    
    print $ofh "data = read.table(\"$counts_matrix_file\", header=T, row.names=1)\n";
    print $ofh "col_ordering = c(" . join(",", @rep_indices) . ")\n";
    print $ofh "rnaseqMatrix = data[,col_ordering]\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)>=5,]\n";
    print $ofh "conditions = factor(c(rep(\"$repA_name\", $num_rep_A), rep(\"$repB_name\", $num_rep_B)))\n";
    print $ofh "\n";
    print $ofh "exp_study = DGEList(counts=rnaseqMatrix, group=conditions)\n";
    print $ofh "exp_study = calcNormFactors(exp_study)\n";
    print $ofh "exp_study = estimateCommonDisp(exp_study)\n";
    print $ofh "exp_study = estimateTagwiseDisp(exp_study)\n";
    print $ofh "et = exactTest(exp_study)\n";
    print $ofh "tTags = topTags(et,n=NULL)\n";
    print $ofh "write.table(tTags[tTags\$table\$PValue <= 0.05,], file=\'$output_prefix.edgeR.output\', sep='\t', quote=F, row.names=T)\n";
        
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
    
    
