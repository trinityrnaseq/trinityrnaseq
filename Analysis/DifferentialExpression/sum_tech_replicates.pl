#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

#######################################################################
#
#  Required:
#
#  --matrix|m <string>          replicates.TMM_normalized.fpkm.matrix
#
#  --samples_file|s <string>    description of sample-to-replicate mapping (tab-delimited)
#
#                             ex.
#                             condA_br1    tech_condA_br1a
#                             condA_br1    tech_condA_br1b
#                             condA_br2    tech_condA_br2a
#                             condA_br2    tech_condA_br2b
#                              
#  --test|t                  write and execute on built-in test data
# 
#
#######################################################################


__EOUSAGE__

    ;


my $matrix_file;
my $use_avg_log_val = 0;
my $samples_file;

my $help_flag;

my $test_flag = 0;

&GetOptions("matrix|m=s" => \$matrix_file,
            "samples_file|s=s" => \$samples_file,

            "help|h" => \$help_flag,

            "test|t" => \$test_flag,

    );

if ($help_flag) {
    die $usage;
}

unless (($matrix_file && $samples_file) || $test_flag) {
    die $usage;
}


main: {


    if ($test_flag) {
        $matrix_file = "test.matrix";
        $samples_file = "test.samples";
        &write_test_data($matrix_file, $samples_file);
    }
    
    

    my $out_prefix = "$matrix_file.sum_tech_reps";
    
    my $out_matrix = $out_prefix . ".matrix";
    my $out_rscript = $out_prefix . ".Rscript";
    
    open (my $ofh, ">$out_rscript") or die "Error, cannot write to $out_rscript";
    print $ofh <<__EORSCRIPT__;

samples = read.table("$samples_file", header=F, check.names=F)
sample_types = as.vector(unique(samples[,1]))
nsamples = length(sample_types)

data = read.table("$matrix_file", header=T, row.names=1, com='', nrows=10000, check.names=F)
classes = sapply(data,class)
data = read.table("$matrix_file", header=T, row.names=1, com='', colClasses=classes, check.names=F)
data = as.matrix(data)
sample_factoring = rep(NA, ncol(data))
sample_expr_matrix = matrix(ncol=nsamples, nrow=nrow(data))
colnames(sample_expr_matrix) = sample_types
rownames(sample_expr_matrix) = rownames(data)

for (i in 1:nsamples) {
    sample_type = sample_types[i]
    rep_indices = samples[,1] \%in% sample_type 
    #print(paste("rep_indices:", rep_indices))
    rep_names = as.vector(samples[rep_indices,2])
    cat("sample_type: ", sample_type, "\n")
    cat("rep_names: ", rep_names, "\n")
    col_indices = colnames(data) \%in% rep_names 
    #print(col_indices);
    if (sum(col_indices) == 0) {
        stop(cat("Error, no columns found matching sample type:", sample_type, "with replicate names:", rep_names))
    }
    if (sum(col_indices) != length(rep_names)) {
        found_colnames = colnames(data)[col_indices]
        missing = rep_names[ ! rep_names \%in% found_colnames]
        cat ("Error, not all replicates accounted for. Found only", sum(col_indices), "columns but have", length(rep_names), "replicates:", rep_names, "for sample:", sample_type, "missing:", missing)
        stop()
    }
    sample_factoring[col_indices] = sample_type
    sum_vals = apply(data[,col_indices, drop=F], 1, sum)
    sample_expr_matrix[,i] = sprintf("%.1f", sum_vals)
}

write.table(sample_expr_matrix, "$out_matrix", quote=F, sep="\t")

__EORSCRIPT__

;

    my $cmd = "R --no-save --no-restore --no-site-file --no-init-file -q < $out_rscript";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }


    print "\n\nDone.  See output: $out_matrix\n\n";
    
    exit(0);
}


####
sub write_test_data {
    my ($test_matrix_file, $test_samples_file) = @_;

    open (my $ofh, ">$test_matrix_file") or die $!;
    print $ofh 
"S1	S2	S3	S4	S5	C1	C2	C3	C4	C5
1	456	726	380	554	184	9	7	3	9	7
2	34	49	43	35	30	18	15	11	23	7
3	22	17	35	17	18	45	127	73	78	93
4	11	6	11	4	7	52	57	44	27	19
5	18	19	12	25	17	193	130	69	137	164
6	41	71	38	59	36	37	45	24	41	16
7	2	11	17	12	10	152	230	137	185	154
8	6	20	26	25	22	2	4	0	1	2
9	12	28	11	13	16	4	15	6	12	12
10	34	31	13	30	25	8	3	10	12	16
";

    close $ofh;


    open ($ofh, ">$test_samples_file") or die $!;
    print $ofh 
"S	S1
S	S2
S	S3
S	S4
S	S5
C	C1
C	C2
C	C3
C	C4
C	C5
";

    return;
}

