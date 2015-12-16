#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use FindBin;
use File::Basename;
use lib ("$FindBin::RealBin/../../PerlLib");
use Data::Dumper;


################################################
# for info on upper quartile normalization, see:
# http://vinaykmittal.blogspot.com/2013/10/fpkmrpkm-normalization-caveat-and-upper.html
################################################

my $usage = <<__EOUSAGE__;



###################################################################################################################
#
#  Required:
#
#  --matrix <string>      FPKM matrix
#
#  --min_expr_in_calc <float>    minimum (exclusive) fpkm value to consider in computing the upper quartile of the distribution.
#                                default: must be > 0
#
####################################################################################################################




__EOUSAGE__


    ;



my $matrix_file;
my $help_flag;
my $min_expr_in_calc = 0;

&GetOptions ( 'h' => \$help_flag,
              'matrix=s' => \$matrix_file,
              'min_expr_in_calc=f' => \$min_expr_in_calc,
              );


if ($help_flag) {
    die $usage;
}


unless ($matrix_file) {
    die $usage;
}



main: {
    
    &upper_quartile_normalize($matrix_file);
    
    exit(0);
}


####
sub upper_quartile_normalize {
    my ($matrix_file) = @_;
    
    my $tmm_norm_script = "__tmp_upper_quart_norm.R";
    open (my $ofh, ">$tmm_norm_script") or die "Error, cannot write to $tmm_norm_script";
    #print $ofh "source(\"$FindBin::RealBin/R/edgeR_funcs.R\")\n";
    
    print $ofh "data = read.table(\"$matrix_file\", header=T, row.names=1, com='')\n";
    print $ofh "get_upper_quartile = function(vec) {\n"
        .      "    vec = vec[vec > $min_expr_in_calc]\n"
        .      "    quantile(vec, 0.75)\n"
        .      "}\n";

    print $ofh "data = as.matrix(data)\n";
    print $ofh "upp_quartiles = apply(data, 2, get_upper_quartile)\n";
    print $ofh "m = sweep(data, MARGIN=2, upp_quartiles, '/')\n";
    print $ofh "mean_upp_quart = mean(upp_quartiles)\n";
    print $ofh "m = m * mean_upp_quart\n";
    
    print $ofh "write.table(m, quote=F, sep=\"\\t\")\n";
    
    close $ofh;
    
    &process_cmd("R --vanilla -q --slave < $tmm_norm_script ");
    
    return;
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret ($ret) ";
    }

    return;
}

