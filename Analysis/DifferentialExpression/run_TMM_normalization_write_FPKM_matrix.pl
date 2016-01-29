#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use FindBin;
use File::Basename;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;
use Data::Dumper;

my $usage = <<__EOUSAGE__;



#################################################################################################
#
#  Required:
#
#  --matrix <string>      matrix of raw read counts (not normalized!)
#  --lengths <string>     extracted from RSEM.(isoforms|genes).results file like so: cat RSEM.isoforms.results | cut -f1,3,4 > lengths.txt
#
#  Optional:
#  --use_eff_length       use 'effective' transcript length instead of effective transcript length.
#                           The effective length represents the number of positions a fragment could be derived from
#                           and is typicaly (total_trans_length - avg_fragment_length + 1)
#                         (for more info, see RSEM paper or cufflinks supp. materials) 
#
#
#  --just_do_TMM_scaling    only rescales the data according to the TMM normalization factors (--lengths param not required)
#
################################################################################################




__EOUSAGE__


    ;



my $matrix_file;
my $help_flag;
my $output_dir;
my $lengths_file;
my $USE_EFF_LENGTH;

my $JUST_DO_TMM_SCALING;

my $dispersion = 0.1;
my $FDR = 0.05;

&GetOptions ( 'h' => \$help_flag,
              'matrix=s' => \$matrix_file,
              'lengths=s' => \$lengths_file,
              'use_eff_length' => \$USE_EFF_LENGTH,
              'just_do_TMM_scaling' => \$JUST_DO_TMM_SCALING,
              );


if ($help_flag) {
    die $usage;
}


unless ($matrix_file && ($lengths_file || $JUST_DO_TMM_SCALING)) { 
    die $usage;
}



main: {

    my %lengths;

    unless ($JUST_DO_TMM_SCALING) {

        %lengths = &get_transcript_lengths($lengths_file, $USE_EFF_LENGTH);
    }
    
    #use Data::Dumper;
    #print Dumper(\%lengths);

    if ($matrix_file !~ /^\//) {
        ## make full path
        $matrix_file = cwd() . "/$matrix_file";
    }
    
    my $tmm_info_file = &run_TMM($matrix_file);
    
    &write_normalized_file($matrix_file, $tmm_info_file, \%lengths);
    
    exit(0);
}


####
sub run_TMM {
    my ($counts_matrix_file) = @_;
    
    my $tmm_norm_script = "__tmp_runTMM.R";
    open (my $ofh, ">$tmm_norm_script") or die "Error, cannot write to $tmm_norm_script";
    #print $ofh "source(\"$FindBin::RealBin/R/edgeR_funcs.R\")\n";
    
    print $ofh "library(edgeR)\n\n";
    
    print $ofh "rnaseqMatrix = read.table(\"$counts_matrix_file\", header=T, row.names=1, com='')\n";
    print $ofh "rnaseqMatrix = round(rnaseqMatrix)\n";
    print $ofh "exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))\n";
    print $ofh "exp_study = calcNormFactors(exp_study)\n";
    
    print $ofh "exp_study\$samples\$eff.lib.size = exp_study\$samples\$lib.size * exp_study\$samples\$norm.factors\n";
    print $ofh "write.table(exp_study\$samples, file=\"$counts_matrix_file.TMM_info.txt\", quote=F, sep=\"\\t\", row.names=F)\n";
    
    close $ofh;

    &process_cmd("R --vanilla -q < $tmm_norm_script");
    
    return("$counts_matrix_file.TMM_info.txt");

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
sub write_normalized_file {
    my ($matrix_file, $tmm_info_file, $seq_lengths_href) = @_;

    my %col_to_eff_lib_size;
    my %col_to_norm_factor;
    open (my $fh, $tmm_info_file) or die "Error, cannot open file $tmm_info_file";
    
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($col, $norm_factor, $eff_lib_size) = ($x[0], $x[2], $x[3]);
        $col_to_eff_lib_size{$col} = $eff_lib_size;
        $col_to_norm_factor{$col} = $norm_factor;
    }
    close $fh;

    my $outfile = ($JUST_DO_TMM_SCALING) ? "$matrix_file.TMM_rescaled"
        : "$matrix_file.TMM_normalized.FPKM";
    
    open (my $ofh, ">$outfile") or die "Error, cannot write to file $outfile";
    
    
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
            print $ofh join("\t", @pos_to_col) . "\n";
        }
        

        my $gene = $x[0];
        my $seq_len;
        unless ($JUST_DO_TMM_SCALING) {
            $seq_len = $seq_lengths_href->{$gene} or die "Error, no seq length for $gene";
        }
        
        print $ofh $gene;
        for (my $i = 1; $i <= $#x; $i++) {
            my $col = $pos_to_col[$i];
            my $adj_col = $col;
            $adj_col =~ s/-/\./g;
                
            my $eff_lib_size = $col_to_eff_lib_size{$col} 
            || $col_to_eff_lib_size{$adj_col}
            || $col_to_eff_lib_size{"X$col"} 
            || die "\n\nError, no eff lib size for [$col]\n\n" . Dumper(\%col_to_eff_lib_size);
            
            my $norm_factor = $col_to_norm_factor{$col}
            || $col_to_norm_factor{$adj_col}
            || $col_to_norm_factor{"X$col"}
            || die "Error, no normalization scaling factor for $col" . Dumper(\%col_to_norm_factor);
            

            my $read_count = $x[$i];
    
            my $converted_val = ($JUST_DO_TMM_SCALING) ? sprintf("%.3f", $read_count * 1/$norm_factor)
                : sprintf("%.2f", $read_count / ($seq_len / 1e3) / ($eff_lib_size / 1e6));
            
            print $ofh "\t$converted_val";
        }
        print $ofh "\n";
    }
    close $ofh;

    return;

}

####
sub get_transcript_lengths {
    my ($file, $use_eff_len) = @_;

    my %lengths;

    open (my $fh, $file) or die "Error, cannot open file $file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        if ($x[2] !~ /^[\d\.]+$/) {
            print STDERR "Skipping line: $_\n";
            next;
        }
        
        
        unless (scalar @x == 3) {
            die "Error, unexpected format of file: $file, should contain format:\n"
                . "accession(tab)length(tab)effective_length\n";
        }
        
        my ($acc, $len, $eff_len) = @x;

        my $len_use = ($use_eff_len) ? $eff_len : $len;
        
        if ( (! $lengths{$acc}) || $len_use > $lengths{$acc}) {
            
            $lengths{$acc} = $len_use;
            
        }
    }
    close $fh;

    return(%lengths);
    
}


