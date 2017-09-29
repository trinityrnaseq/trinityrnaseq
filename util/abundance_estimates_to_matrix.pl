#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Carp;

my $usage = <<__EOUSAGE__;

####################################################################################
#
# Usage:  $0 --est_method <method>  sample1.results sample2.results ...
#
#      or  $0 --est_method <method> --quant_files file.listing_target_files.txt
#
#      Note, if only a single input file is given, it's expected to contain the paths to all the target abundance estimation files.
#
# Required:
#            
#  --est_method <string>           RSEM|eXpress|kallisto|salmon  (needs to know what format to expect)
#
#  --gene_trans_map <string>           the gene-to-transcript mapping file. (if you don't want gene estimates, indicate 'none'.
#
#
# Options:
#
#  --cross_sample_norm <string>         TMM|UpperQuartile|none   (default: TMM)
#
#  --name_sample_by_basedir             name sample column by dirname instead of filename
#      --basedir_index <int>            default(-2)
#
#  --out_prefix <string>                default: value for --est_method
#
#  --quant_files <string>              file containing a list of all the target files.
#
######################################################################################


__EOUSAGE__

    ;


my $help_flag;
my $est_method;
my $val_type;
my $cross_sample_norm = "TMM";
my $name_sample_by_basedir = 0;
my $out_prefix;
my $basedir_index = -2;
my $quant_files = "";
my $gene_trans_map_file;

&GetOptions('help|h' => \$help_flag,
            'est_method=s' => \$est_method,
            
            'cross_sample_norm=s' => \$cross_sample_norm,
            'name_sample_by_basedir' => \$name_sample_by_basedir,
            'out_prefix=s' => \$out_prefix,
            'basedir_index=i' => \$basedir_index,
            
            'quant_files=s' => \$quant_files,
            'gene_trans_map=s' => \$gene_trans_map_file,
    );



if ($help_flag) { die $usage; }

unless ($est_method && (@ARGV || $quant_files)) {
    die $usage;
}

unless ($gene_trans_map_file) {
    die "Error, specify gene-to-trans map file via: --gene_trans_map, or indicate 'none' if you dont want gene estimates";
}

unless ($est_method =~ /^(RSEM|eXpress|kallisto|salmon)/i) {
    die "Error, dont recognize --est_method $est_method ";
}
unless ($cross_sample_norm =~ /^(TMM|UpperQuartile|none)$/i) {
    die "Error, dont recognize --cross_sample_norm $cross_sample_norm ";
}

unless ($out_prefix) {
    $out_prefix = $est_method;
}

my @files;

if ($quant_files) {
    # allow for a file listing the various files.
    @files = `cat $quant_files`;
    chomp @files;
}
elsif (@ARGV) {
    @files = @ARGV;
}
else {
    die $usage;
}



=data_formats

## RSEM:

0       transcript_id
1       gene_id
2       length
3       effective_length
4       expected_count
5       TPM
6       FPKM
7       IsoPct


## eXpress v1.5:

1       target_id
2       length
3       eff_length
4       tot_counts
5       uniq_counts
6       est_counts
7       eff_counts
8       ambig_distr_alpha
9       ambig_distr_beta
10      fpkm
11      fpkm_conf_low
12      fpkm_conf_high
13      solvable
14      tpm


## kallisto:
0       target_id
1       length
2       eff_length
3       est_counts
4       tpm


## salmon:
0       Name
1       Length
2       EffectiveLength
3       TPM
4       NumReads

=cut
    
    ;

my ($acc_field, $counts_field, $fpkm_field, $tpm_field);

if ($est_method =~ /^rsem$/i) {
    $acc_field = 0;
    $counts_field = "expected_count";
    $fpkm_field = "FPKM";
    $tpm_field = "TPM";
}
elsif ($est_method =~ /^express$/i) {  # as of v1.5
    $acc_field = "target_id";
    $counts_field = "eff_counts";
    $fpkm_field = "fpkm";
    $tpm_field = "tpm";
}
elsif ($est_method =~ /^kallisto$/i) {
    $acc_field = "target_id";
    $counts_field = "est_counts";
    $fpkm_field = "tpm";
    $tpm_field = "tpm";
}
elsif ($est_method =~ /^salmon/) {
    $acc_field = "Name";
    $counts_field = "NumReads";
    $fpkm_field = "TPM";
    $tpm_field = "TPM";
}
else {
    die "Error, dont recognize --est_method [$est_method] ";
}

main: {
    
    my %data;

    my %sum_sample_counts;
    
    foreach my $file (@files) {
        print STDERR "-reading file: $file\n";
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>; 
        chomp $header;
        my %fields = &parse_field_positions($header);
        #use Data::Dumper; print STDERR Dumper(\%fields);
        while (<$fh>) {
            chomp;
            
            my @x = split(/\t/);
            my $acc = $x[ $fields{$acc_field} ];
            my $count = $x[ $fields{$counts_field} ];
            my $fpkm = $x[ $fields{$fpkm_field} ];
            my $tpm = $x[ $fields{$tpm_field} ];

            $data{$acc}->{$file}->{count} = $count;
            $data{$acc}->{$file}->{FPKM} = $fpkm;
            $data{$acc}->{$file}->{TPM} = $tpm;
            
            # capture sample total counts
            $sum_sample_counts{$file} += $count;
        }
        close $fh;
    }

    my %column_header_to_filename;
    my @filenames = @files;
    foreach my $file (@filenames) {
        my $column_header;
        if ($name_sample_by_basedir) {
            my @path = split(m|/|, $file);
            $column_header = $path[$basedir_index];
        }
        else {
            $column_header = basename($file);
            
        }
        $column_header_to_filename{$column_header} = $file;
        $file = $column_header; # update the @filenames
    }
    print STDERR "\n\n* Outputting combined matrix.\n\n";
    
    my $counts_matrix_file = "$out_prefix.isoform.counts.matrix";
    my $TPM_matrix_file = "$out_prefix.isoform.TPM.not_cross_norm";
    open (my $ofh_counts, ">$counts_matrix_file") or die "Error, cannot write file $counts_matrix_file";
    open (my $ofh_TPM, ">$TPM_matrix_file") or die "Error, cannot write file $TPM_matrix_file";
    
    { # check to see if they're unique
        my %filename_map = map { + $_ => 1 } @filenames;
        if (scalar keys %filename_map != scalar @filenames) {
            die "Error, the column headings: @filenames are not unique.  Should you consider using the --name_sample_by_basedir parameter?";
        }
    }
    

    # clean up matrix headers
    foreach my $file (@filenames) {
        # also, get rid of the part of the filename that RSEM adds
        $file =~ s/\.(genes|isoforms)\.results$//;
    }
    

    print $ofh_counts join("\t", "", @filenames) . "\n";
    print $ofh_TPM join("\t", "", @filenames) . "\n";

    foreach my $acc (keys %data) {
        
        print $ofh_counts "$acc";
        print $ofh_TPM "$acc";
        
        foreach my $file (@files) {

            my $count = $data{$acc}->{$file}->{count};
            unless (defined $count) {
                $count = "NA";
            }
            my $tpm = $data{$acc}->{$file}->{TPM};
            if (defined $tpm) {
                $tpm = $tpm/1;
            }
            else {
                $tpm = "NA";
            }

            print $ofh_counts "\t$count";
            print $ofh_TPM "\t$tpm";
        }
        
        print $ofh_counts "\n";
        print $ofh_TPM "\n";

    }
    close $ofh_counts;
    close $ofh_TPM;

    ## process gene counts as per txImport-style (Soneson et al. F1000, 2016)
    my $gene_counts_file = "$out_prefix.gene.counts.matrix";
    my $gene_tpm_file = "$out_prefix.gene.TPM.not_cross_norm";
    
    if ($gene_trans_map_file) {
        my %gene_to_trans;
        {
            open(my $fh, $gene_trans_map_file) or die "Error, cannot open file $gene_trans_map_file";
            while (<$fh>) {
                chomp;
                my ($gene, $trans) = split(/\s+/);
                push (@{$gene_to_trans{$gene}}, $trans);
            }
            close $fh;
        }

        open(my $ofh_genecounts, ">$gene_counts_file") or die "Error, cannot write to $gene_counts_file";
        print $ofh_genecounts "\t" . join("\t", @filenames) . "\n";
        
        open(my $ofh_genetpm, ">$gene_tpm_file") or die "Error, cannot write to $gene_tpm_file";
        print $ofh_genetpm "\t" . join("\t", @filenames) . "\n";
        
        foreach my $gene (sort keys %gene_to_trans) {
            my @tpm_vals = ($gene);
            my @count_vals = ($gene);
            foreach my $file (@filenames) {
                my $gene_tpm = 0;
                # sum up gene tpm from isoform tpms
                foreach my $trans (@{$gene_to_trans{$gene}}) {
                    my $trans_tpm = $data{$trans}->{ $column_header_to_filename{$file} }->{TPM};
                    unless (defined $trans_tpm) {
                        confess "Error, no TPM value specified for transcript [$trans] of gene [$gene] for sample $file";
                    }
                    $gene_tpm += $trans_tpm;
                }
                push (@tpm_vals, $gene_tpm);
                my $gene_count = $gene_tpm / 1e6 * $sum_sample_counts{ $column_header_to_filename{$file} };
                push (@count_vals, $gene_count);
            }
            
            print $ofh_genetpm join("\t", @tpm_vals) . "\n";
            print $ofh_genecounts join("\t", @count_vals);
        }
        close $ofh_genetpm;
        close $ofh_genecounts;
    }
    if (scalar @files > 1) {
        ## more than one sample 
        
        &perform_cross_sample_norm($TPM_matrix_file, "$out_prefix.isoform");
        if ($gene_trans_map_file) {
            &perform_cross_sample_norm($gene_tpm_file, "$out_prefix.gene");
        }
        
    }
    else {
        unless (scalar @files == 1) { 
            die "Error, no target samples. Shouldn't get here.";
        }
        print STDERR "Warning, only one sample, so not performing cross-sample normalization\n";
        print STDERR "Done.\n\n";
    }
    
    exit(0);
}


####
sub perform_cross_sample_norm {
    my ($tpm_matrix_file, $out_prefix_name) = @_;
    
    if ($cross_sample_norm =~ /^TMM$/i) {
        my $cmd = "$FindBin::RealBin/support_scripts/run_TMM_scale_matrix.pl --matrix $tpm_matrix_file > $out_prefix_name.$cross_sample_norm.EXPR.matrix";
        &process_cmd($cmd);
    }
    elsif ($cross_sample_norm =~ /^UpperQuartile$/) {
        my $cmd = "$FindBin::RealBin/support_scripts/run_UpperQuartileNormalization_matrix.pl --matrix $tpm_matrix_file > $out_prefix_name.$cross_sample_norm.EXPR.matrix";
        &process_cmd($cmd);
    }
    elsif ($cross_sample_norm =~ /^none$/i) {
        print STDERR "-not performing cross-sample normalization.\n";
    }

}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR $cmd;
    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}
        

####
sub parse_field_positions {
    my ($header) = @_;

    my %field_pos;
    my @fields = split(/\s+/, $header);
    for (my $i = 0; $i <= $#fields; $i++) {
        $field_pos{$fields[$i]} = $i;
        $field_pos{$i} = $i; # for fixed column assignment
    }

    return(%field_pos);
}
