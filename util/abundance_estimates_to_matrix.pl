#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --est_method <method>  sample1.results sample2.results ...
# Required:
#
#  --est_method <string>           RSEM|eXpress  (needs to know what format to expect)
#
#
# Options:
#
#  --cross_sample_fpkm_norm <string>    TMM|UpperQuartile|none   (default: TMM)
#
#  --name_sample_by_basedir             name sample column by dirname instead of filename
#      --basedir_index <int>            default(-2)
#
#  --out_prefix <string>                default: 'matrix'
#
############################################################


__EOUSAGE__

    ;


my $help_flag;
my $est_method;
my $cross_sample_fpkm_norm = "TMM";
my $name_sample_by_basedir = 0;
my $out_prefix = "matrix";
my $basedir_index = -2;

&GetOptions('help|h' => \$help_flag,
            'est_method=s' => \$est_method,
            
            'cross_sample_fpkm_norm=s' => \$cross_sample_fpkm_norm,
            'name_sample_by_basedir' => \$name_sample_by_basedir,
            'out_prefix=s' => \$out_prefix,
            
            'basedir_index=i' => \$basedir_index,
            
            );


unless ($est_method && @ARGV) {
    die $usage;
}

unless ($est_method =~ /^(RSEM|eXpress)$/i) {
    die "Error, dont recognize --est_method $est_method ";
}
unless ($cross_sample_fpkm_norm =~ /^(TMM|UpperQuartile|none)$/i) {
    die "Error, dont recognize --cross_sample_fpkm_norm $cross_sample_fpkm_norm ";
}

my @files = @ARGV;

if (scalar @files == 1) {

    if (-s $files[0]) {
        # allow for a file listing the various files.
        @files = `cat $files[0]`;
        chomp @files;
    }
    else {
        die $usage;
    }
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


## eXpress:

0       bundle_id
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

=cut
    
    ;

my ($acc_field, $counts_field, $fpkm_field);

if ($est_method =~ /^rsem$/i) {
    $acc_field = 0;
    $counts_field = 4;
    $fpkm_field = 6;
}
elsif ($est_method =~ /^express$/i) {
    $acc_field = 1;
    $counts_field = 7;
    $fpkm_field = 10;
}
else {
    die "Error, dont recognize --est_method $est_method ";
}

main: {
    
    my %data;
    
    foreach my $file (@files) {
        print STDERR "-reading file: $file\n";
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>; # ignore it
        while (<$fh>) {
            chomp;
            
            my @x = split(/\t/);
            my $acc = $x[$acc_field];
            my $count = $x[$counts_field];
            my $fpkm = $x[$fpkm_field];
            
            $data{$acc}->{$file}->{count} = $count;
            $data{$acc}->{$file}->{fpkm} = $fpkm;
        }
        close $fh;
    }
    
    my @filenames = @files;
    foreach my $file (@filenames) {
        if ($name_sample_by_basedir) {
            my @path = split(m|/|, $file);
            $file = $path[$basedir_index];
        }
        else {
            $file = basename($file);
        }
    }
    print STDERR "\n\n* Outputting combined matrix.\n\n";
    
    my $counts_matrix_file = "$out_prefix.counts.matrix";
    my $fpkm_matrix_file = "$out_prefix.not_cross_norm.fpkm.tmp";
    open (my $ofh_counts, ">$counts_matrix_file") or die "Error, cannot write file $counts_matrix_file";
    open (my $ofh_fpkm, ">$fpkm_matrix_file") or die "Error, cannot write file $fpkm_matrix_file";

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
    print $ofh_fpkm join("\t", "", @filenames) . "\n";

    foreach my $acc (keys %data) {
        
        print $ofh_counts "$acc";
        print $ofh_fpkm "$acc";
        
        foreach my $file (@files) {

            my $count = $data{$acc}->{$file}->{count};
            unless (defined $count) {
                $count = "NA";
            }
            my $fpkm = $data{$acc}->{$file}->{fpkm};
            unless (defined $fpkm) {
                $fpkm = "NA";
            }

            print $ofh_counts "\t$count";
            print $ofh_fpkm "\t$fpkm";
        }
        
        print $ofh_counts "\n";
        print $ofh_fpkm "\n";

    }
    close $ofh_counts;
    close $ofh_fpkm;
    
    if ($cross_sample_fpkm_norm =~ /^TMM$/i) {
        my $cmd = "$FindBin::Bin/support_scripts/run_TMM_scale_matrix.pl --matrix $fpkm_matrix_file > $out_prefix.$cross_sample_fpkm_norm.fpkm.matrix";
        &process_cmd($cmd);
    }
    elsif ($cross_sample_fpkm_norm =~ /^UpperQuartile$/) {
        my $cmd = "$FindBin::Bin/support_scripts/run_UpperQuartileNormalization_matrix.pl --matrix $fpkm_matrix_file > $out_prefix.$cross_sample_fpkm_norm.fpkm.matrix";
        &process_cmd($cmd);
    }
    elsif ($cross_sample_fpkm_norm =~ /^none$/i) {
        print STDERR "-not performing cross-sample normalization.\n";
    }
    
    print STDERR "Done.\n\n";
    
    exit(0);
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
        
