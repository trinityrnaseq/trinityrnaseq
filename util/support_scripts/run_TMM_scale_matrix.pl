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

my $usage = <<__EOUSAGE__;



#################################################################################################
#
#  Required:
#
#  --matrix <string>      matrix of raw read counts (not normalized!)
#
################################################################################################




__EOUSAGE__


    ;



my $matrix_file;
my $help_flag;

&GetOptions ( 'h' => \$help_flag,
              'matrix=s' => \$matrix_file,
              );


if ($help_flag) {
    die $usage;
}


unless ($matrix_file) {
    die $usage;
}



main: {
    
    my $tmm_info_file = &run_TMM($matrix_file);
    
    &write_normalized_file($matrix_file, $tmm_info_file);
    
    exit(0);
}


####
sub run_TMM {
    my ($counts_matrix_file) = @_;
    
    my $tmm_norm_script = "__tmp_runTMM.R";
    open (my $ofh, ">$tmm_norm_script") or die "Error, cannot write to $tmm_norm_script";
    #print $ofh "source(\"$FindBin::RealBin/R/edgeR_funcs.R\")\n";
    
    print $ofh "library(edgeR)\n\n";
    
    print $ofh "rnaseqMatrix = read.table(\"$counts_matrix_file\", header=T, row.names=1, com='', check.names=F)\n";
    print $ofh "rnaseqMatrix = as.matrix(rnaseqMatrix)\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))\n";
    print $ofh "exp_study = calcNormFactors(exp_study)\n";
    
    print $ofh "exp_study\$samples\$eff.lib.size = exp_study\$samples\$lib.size * exp_study\$samples\$norm.factors\n";
    print $ofh "write.table(exp_study\$samples, file=\"$counts_matrix_file.TMM_info.txt\", quote=F, sep=\"\\t\", row.names=F)\n";
    
    close $ofh;

    &process_cmd("R --vanilla -q < $tmm_norm_script 1>&2 ");
    
    my $tmm_matrix = "$counts_matrix_file.TMM_info.txt";

    unless (-s $tmm_matrix) {
        confess "Error, TMM matrix $tmm_matrix was not generated.  Be sure edgeR is installed and see additional error messages above for other helpful info";
    }
    
    return($tmm_matrix);
    
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

####
sub write_normalized_file {
    my ($matrix_file, $tmm_info_file) = @_;

    my %col_to_eff_lib_size;
    my %col_to_norm_factor;
    open (my $fh, $tmm_info_file) or die "Error, cannot open file $tmm_info_file";
    
    my %bloom_to_col;

    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($col, $norm_factor, $eff_lib_size) = ($x[0], $x[2], $x[3]);
        
        $col =~ s/\"//g;
        
        my $bloom = $col;
        $bloom =~ s/\W/$;/g;
        
        $col_to_eff_lib_size{$col} = $eff_lib_size;
        $col_to_norm_factor{$col} = $norm_factor;
    
        if ($bloom ne $col) {

            if (exists $bloom_to_col{$bloom}) {
                die "Error, already stored $bloom_to_col{$bloom} for $bloom, but trying to also store $col here... Ensure column names are unique according to non-word characters.";
            }

            $col_to_eff_lib_size{$bloom} = $eff_lib_size;
            $col_to_norm_factor{$bloom} = $norm_factor;
        }
        
    }
    close $fh;
    
    open ($fh, $matrix_file);
    $header = <$fh>;
    chomp $header;
    my @pos_to_col = split(/\t/, $header);
    my $check_column_ordering_flag = 0;
    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        unless ($check_column_ordering_flag) {
            if (scalar(@x) == scalar(@pos_to_col) + 1) {
                ## header is offset, as is acceptable by R
                ## not acceptable here.  fix it:
                unshift (@pos_to_col, "");
            }
            $check_column_ordering_flag = 1;
            print join("\t", @pos_to_col) . "\n";
        }
        

        my $gene = $x[0];
        
        print $gene;
        for (my $i = 1; $i <= $#x; $i++) {
            my $col = $pos_to_col[$i];
            
            $col =~ s/\"//g;
            
            my $adj_col = $col;
            $adj_col =~ s/-/\./g;
            
            my $bloom = $col;
            $bloom =~ s/\W/$;/g;

            my $eff_lib_size = $col_to_eff_lib_size{$col} 
            || $col_to_eff_lib_size{$bloom}
            || $col_to_eff_lib_size{$adj_col}
            || $col_to_eff_lib_size{"X$col"} 
            || die "Error, no eff lib size for [$col] or [$bloom] or [$adj_col] or [\"X$col\"]" . Dumper(\%col_to_eff_lib_size);
            
            my $norm_factor = $col_to_norm_factor{$col}
            || $col_to_norm_factor{$bloom}
            || $col_to_norm_factor{$adj_col}
            || $col_to_norm_factor{"X$col"}
            || die "Error, no normalization scaling factor for $col" . Dumper(\%col_to_norm_factor);
            

            my $read_count = $x[$i];
    
            my $converted_val = sprintf("%.3f", $read_count * 1/$norm_factor);
            
            print "\t$converted_val";
        }
        print "\n";
    }
    
    return;

}
