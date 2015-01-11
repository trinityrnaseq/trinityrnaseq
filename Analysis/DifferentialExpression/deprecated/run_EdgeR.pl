#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use FindBin;
use File::Basename;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $usage = <<__EOUSAGE__;



#################################################################################################
#
#  Required:
#
#  --matrix <string>             matrix of raw read counts (not normalized!)
#  --transcript_lengths <string>   extracted from RSEM.isoforms.results file like so: cat RSEM.isoforms.results | cut -f1,3,4 > lengths.txt
#
#  Optional:
#
#  --reference <string>      specify a reference column to which the others should be compared.
#  
#  --output                  name of directory to place outputs (default: edgeR.\$pid.dir)
#
#  --just_TMM                just do the TMM normalization, don't run the edgeR comparisons.
#
#  --dispersion              edgeR dispersion value (default: 0.1)   set to 0 for poisson (sometimes breaks...)
#
#  --FDR                     false discovery rate (default: 0.05)
#
#  --no_eff_length           use total transcript length instead of effective transcript length.
#
################################################################################################




__EOUSAGE__


    ;



my $matrix_file;
my $reference;
my $help_flag;
my $output_dir;
my $JUST_TMM;
my $transcript_lengths_file;
my $NO_EFF_LENGTH;

my $dispersion = 0.1;
my $FDR = 0.05;

&GetOptions ( 'h' => \$help_flag,
              'matrix=s' => \$matrix_file,
              'transcript_lengths=s' => \$transcript_lengths_file,
              'reference=s' => \$reference,
              'output=s' => \$output_dir,
              'just_TMM' => \$JUST_TMM,
              'dispersion=f' => \$dispersion,
              'FDR=f' => \$FDR,
              'no_eff_length' => \$NO_EFF_LENGTH,
              );


if ($help_flag) {
    die $usage;
}


unless ($matrix_file && $transcript_lengths_file) { 
    die $usage;
}



main: {

    my %trans_lengths = &get_transcript_lengths($transcript_lengths_file, $NO_EFF_LENGTH);
    
    #use Data::Dumper;
    #print Dumper(\%trans_lengths);

    if ($matrix_file !~ /^\//) {
        ## make full path
        $matrix_file = cwd() . "/$matrix_file";
    }
    
    my $edgeR_dir = $output_dir || "edgeR.$$.dir";
    mkdir($edgeR_dir) or die "Error, cannot mkdir $edgeR_dir";
    chdir $edgeR_dir or die "Error, cannot cd to $edgeR_dir";
    
    

    my %filenames = &prepare_edgeR_inputs($matrix_file);
    my @cols = keys %filenames;
    
    my $tmm_info_file = &run_TMM(\%filenames);  ## compute effective library sizes, not used directly here, but will be used in other apps.
    
    my $fpkm_matrix_file = &write_normalized_fpkm_file($matrix_file, $tmm_info_file, \%trans_lengths);
    
    if ($JUST_TMM) {
        print STDERR "-only running TMM, not identifying diff. expressed genes. Stopping now.\n";
        exit(0);
    }
        
   

    my @comparisons;
    
    if ($reference) {
        
        my @comparators = grep { $_ ne $reference } @cols;

        foreach my $comparator (@comparators) {
            my $result_file = &run_edgeR($reference, $comparator, \%filenames);
            push (@comparisons, [$reference, $comparator, $result_file]);
        }
    }
    else {
        ## running all-vs-all
        
        for (my $i = 0; $i < $#cols; $i++) {


            for (my $j = $i + 1; $j <= $#cols; $j++) {

                
                my $result_file = &run_edgeR($cols[$i], $cols[$j], \%filenames);
                push (@comparisons, [$cols[$i], $cols[$j], $result_file]);
            }
        }
    }

    
    ## Summarize all results from the pairwise comparisons in a single output file.

    open (my $ofh, ">all_diff_expression_results.txt") or die $!;
    
    foreach my $comparison (@comparisons) {
        my ($sampleA, $sampleB, $result_file) = @$comparison;
        open (my $fh, $result_file) or die $!;
        my $header = <$fh>;
        print $ofh "#sampleA\tsampleB\ttranscript\t$header";
        while (<$fh>) {
            print $ofh join("\t", $sampleA, $sampleB, $_);
        }
        close $fh;
    }
    close $ofh;

    print STDERR "\n\nAll diff expressed results written to: all_diff_expression_results.txt\n\n";
            
    ## Generate diff expression report summary
    

    exit(0);
}


####
sub prepare_edgeR_inputs {
    my ($matrix_file) = @_;
    
    open (my $fh, $matrix_file) or die $!;
    my $col_header = <$fh>;
    chomp $col_header;
    
    my @cols = split(/\t/, $col_header);
    shift @cols; # rid gene name
    
    my %fhs;
    my %filenames;
    foreach my $col (@cols) {
        my $colname = $col;
        $colname =~ s/\W/_/g;
        my $dat_file = "$colname.dat";
        open (my $ofh, ">$dat_file") or die "Error, cannot write to $dat_file";
        $fhs{$col} = $ofh;
        $filenames{$col} = $dat_file;
    }
    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $gene_name = shift @x;
        foreach my $col (@cols) {
            my $ofh = $fhs{$col};
            my $val = shift @x;
            print $ofh "$gene_name\t$val\n";
        }
    }
    close $fh;
    
    foreach my $ofh (values %fhs) {
        close $ofh;
    }
    
    return(%filenames);
}



####
sub run_TMM {
    my ($filenames_href) = @_;
        
    my $data_file_list = "all_data_files.list";
    open (my $ofh, ">$data_file_list") or die $!;
    print $ofh "files\tgroup\tdescription\n";
    
    foreach my $col (keys %$filenames_href) {
        my $filename = $filenames_href->{$col};
        
        print $ofh join("\t", $filename, $col, "$col sample") . "\n";
    }
    close $ofh;

    my $tmm_norm_script = "__tmp_runTMM.R";
    open ($ofh, ">$tmm_norm_script") or die "Error, cannot write to $tmm_norm_script";
    print $ofh "source(\"$FindBin::Bin/R/edgeR_funcs.R\")\n";
    print $ofh "myDGEList = target_list_file_to_DGEList(\"$data_file_list\")\n";
    print $ofh "myDGEList\$samples\$eff.lib.size = myDGEList\$samples\$lib.size * myDGEList\$samples\$norm.factors\n";
    print $ofh "write.table(myDGEList\$samples, file=\"TMM_info.txt\", quote=F, sep=\"\\t\", row.names=F)\n";
    
    close $ofh;

    &process_cmd("R --vanilla -q < $tmm_norm_script");
    
    return("TMM_info.txt");

}



####
sub run_edgeR {
    my ($col_A, $col_B, $filenames_href) = @_;
    
    my $R_script = "__tmp_run_edgeR.$$.R";
    my $target_file = "__tmp_targets.$$.dat";
        
    my $compare_name = "${col_A}_vs_${col_B}";
    
    open (my $ofh, ">$target_file") or die "Error, cannot write to $target_file";
    print $ofh "files\tgroup\tdescription\n";
    print $ofh join("\t", $filenames_href->{$col_A}, $col_A, "$col_A sample") . "\n";
    print $ofh join("\t", $filenames_href->{$col_B}, $col_B, "$col_B sample") . "\n";
    close $ofh;
    
    open ($ofh, ">$R_script") or die "Error, cannot write $R_script";
    print $ofh "source(\"$FindBin::Bin/R/edgeR_funcs.R\")\n";
    print $ofh "myDGEList = target_list_file_to_DGEList(\"$target_file\")\n";
    print $ofh "postscript(\"$compare_name.eps\")\n";
    print $ofh "top_entries_table = edgeR_DE_analysis(myDGEList, dispersion=$dispersion, FDR=$FDR); # set dispersion to zero for Poisson-equivalent results.\n";
    print $ofh "dev.off()\n";
    print $ofh "write.table(top_entries_table, quote=F, file=\"$compare_name.results.txt\", sep=\"\\t\")\n";
    
    close $ofh;
    
    my $cmd = "R --vanilla -q < $R_script";
    &process_cmd($cmd);
    

    return ("$compare_name.results.txt");
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
sub write_normalized_fpkm_file {
    my ($matrix_file, $tmm_info_file, $seq_lengths_href) = @_;

    my %col_to_eff_lib_size;
    open (my $fh, $tmm_info_file) or die "Error, cannot open file $tmm_info_file";
    
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($col, $eff_lib_size) = ($x[1], $x[5]);
        $col_to_eff_lib_size{$col} = $eff_lib_size;
    }
    close $fh;

    my $normalized_fpkm_file = "matrix.TMM_normalized.FPKM";
    open (my $ofh, ">$normalized_fpkm_file") or die "Error, cannot write to file $normalized_fpkm_file";
    
    
    open ($fh, $matrix_file);
    $header = <$fh>;
    print $ofh $header;
    chomp $header;
    my @pos_to_col = split(/\t/, $header);
    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $gene = $x[0];
        my $seq_len = $seq_lengths_href->{$gene} or die "Error, no seq length for $gene";
        
        print $ofh $gene;
        for (my $i = 1; $i <= $#x; $i++) {
            my $col = $pos_to_col[$i];
            my $eff_lib_size = $col_to_eff_lib_size{$col} or die "Error, no eff lib size for $col";
            
            my $read_count = $x[$i];
            
            #print STDERR "gene: $gene, read_count: $read_count, seq_len: $seq_len, eff_lib_size: $eff_lib_size\n";
            my $fpkm = ($seq_len > 0) 
                ? $read_count / ($seq_len / 1e3) / ($eff_lib_size / 1e6) 
                : 0;
            
            $fpkm = sprintf("%.2f", $fpkm);

            print $ofh "\t$fpkm";
        }
        print $ofh "\n";
    }
    close $ofh;

    return($normalized_fpkm_file);

}

####
sub get_transcript_lengths {
    my ($file, $no_eff_len) = @_;

    my %trans_lengths;

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

        my $len_use = ($no_eff_len) ? $len : $eff_len;
        
        $trans_lengths{$acc} = $len_use;
        
    }
    close $fh;

    return(%trans_lengths);
    
}


