#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use FindBin;
use Cwd;

######################################################
## Set to base directory of the Trinity installation:
my $BASEDIR = "$FindBin::RealBin/../";
######################################################

my $usage = <<__EOUSAGE__;


##########################################################################################################
#
#  Required:
#
#  --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                                   # note, include Trinity parameters in the sample file like so:
#                                   --CPU=6
#                                   --max_memory=10G
#                                   --seqType=fq
#                                   --SS_lib_type=RF
#
###########################################################################################################

__EOUSAGE__

    ;


my $help_flag;
my $read_samples_descr_file;

&GetOptions ( 'h' => \$help_flag,
              'samples_file=s' => \$read_samples_descr_file,


);

if ($help_flag || ! $read_samples_descr_file) {
    die $usage;
}

my $workdir = cwd();

my @PARAMS;
my %conditions_to_read_info = &parse_sample_descriptions($read_samples_descr_file, \@PARAMS);

my @conditions = sort keys %conditions_to_read_info;


my @left_fq_files;
my @right_fq_files;

foreach my $condition (@conditions) {
    
    my $replicates_href = $conditions_to_read_info{$condition};
    
    my @replicates = keys %$replicates_href;
    foreach my $replicate (@replicates) {
        my ($left_fq_file, $right_fq_file) = @{$replicates_href->{$replicate}};

        if (-s $left_fq_file) {
            push (@left_fq_files, $left_fq_file);
        }
        else {
            die "Error, cannot locate file: [$left_fq_file] ";
        }

        if ($right_fq_file) {
            if (-s $right_fq_file) {
                push (@right_fq_files, $right_fq_file);
            }
            else {
                die "Error, cannot locate file: [$right_fq_file]";
            }
        }
    }
}



## Run Trinity:
my $cmd = "";

if (@left_fq_files && @right_fq_files) {

    $cmd = "$BASEDIR/Trinity --left " . join(",", @left_fq_files) . " --right " . join(",", @right_fq_files);
}
else {
    # run left as single
    $cmd = "$BASEDIR/Trinity --single " . join(",", @left_fq_files);
}

$cmd .= " " . join(" ", @PARAMS);

print STDERR "CMD: $cmd\n\n";

my $ret = system($cmd);

if ($ret) {
    die "Error, CMD: $cmd died with ret: $ret";
}


exit(0);


####
sub parse_sample_descriptions {
    my ($read_samples_descr_file, $params_aref) = @_;

    my %samples_descr;
    
    
    open (my $fh, $read_samples_descr_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        s/^\s+|\s+$//g;
        chomp;
        if (/^\-/) {
            ## a parameter
            push (@$params_aref, $_);
        }
        else {
            my @x = split(/\t/);
            
            my ($condition, $replicate, $reads_left, $reads_right) = @x;
            
            $samples_descr{$condition}->{$replicate} = [$reads_left, $reads_right];
        }
    }
    
    close $fh;
    
    return(%samples_descr);
    
    
}
