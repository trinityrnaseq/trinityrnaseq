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
#  --DE_method <string>            DE method (options: edgeR, DESeq2)
#
#  --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                                   # note, Trinity-specific parameter settings should be included in the samples_file like so:
#                                   # (only --max_memory is absolutely required, since defaults exist for the other settings)
#                                   --CPU=6
#                                   --max_memory=10G
#                                   --seqType=fq
#                                   --SS_lib_type=RF
#
#
#
#
#  Optional:
#
#  -I                              Interactive mode, waits between commands.
#
###########################################################################################################

__EOUSAGE__

    ;


my $help_flag;
my $read_samples_descr_file;
my $PAUSE_EACH_STEP = 0;
my $DE_method;

&GetOptions ( 'h' => \$help_flag,
              'samples_file=s' => \$read_samples_descr_file,
              'I' => \$PAUSE_EACH_STEP,
              'DE_method=s' => \$DE_method,

);

if ($help_flag) {
    die $usage;
}

unless ($read_samples_descr_file && $DE_method) {
    die $usage;
}

unless ($read_samples_descr_file =~ /^\//) {
    $read_samples_descr_file = cwd() . "/$read_samples_descr_file";
}


{
    ## Check for required software
    my @needed_tools = qw (R bowtie bowtie-build); # Trinity.pl, RSEM, and samtools are set by relative paths.
    my $missing_flag = 0;
    foreach my $prog (@needed_tools) {
        my $path = `which $prog`;
        unless ($path =~ /\w/) {
            print STDERR "\n** ERROR, cannot find path to required software: $prog **\n";
            $missing_flag = 1;
        }
    }
    if ($missing_flag) {
        die "\nError, at least one required software tool could not be found. Please install tools and/or adjust your PATH settings before retrying.\n";
    }
}


my $workdir = cwd();

my %PARAMS; 
my %conditions_to_read_info = &parse_sample_descriptions($read_samples_descr_file, \%PARAMS);

my @conditions = sort keys %conditions_to_read_info;


##################
## run DE analysis
##################

foreach my $target_type ("trans", "genes") {
 
    my $edgeR_dir = "${DE_method}_${target_type}";

    my $cmd = "$BASEDIR/Analysis/DifferentialExpression/run_DE_analysis.pl "
        . " --matrix Trinity_${target_type}.counts.matrix "
        . " --method $DE_method "
        . " --samples_file $read_samples_descr_file "
        . " --output $edgeR_dir ";
    
    &process_cmd($cmd, "Running edgeR for $target_type") unless (-d $edgeR_dir);

    
    chdir $edgeR_dir or die "Error, cannot cd to $edgeR_dir";
    
    ## extract the diff. expressed transcripts.
    $cmd = "$BASEDIR/Analysis/DifferentialExpression/analyze_diff_expr.pl "
        . " --matrix ../Trinity_${target_type}.TMM.fpkm.matrix --samples $read_samples_descr_file ";
    

    if (exists $PARAMS{"-P"}) {
        $cmd .= " -P " . $PARAMS{"-P"};
    }
    if (exists $PARAMS{"-C"}) {
        $cmd .= " -C " . $PARAMS{"-C"};
    }
    
    &process_cmd($cmd, "Running analysis of DE $target_type, hierarchically clustering transcripts and samples");
    
    chdir $workdir or die "Error, cannot cd back to $workdir";
}


print "Done.\n";
    



exit(0);


####
sub process_cmd {
    my ($cmd, $msg) = @_;


    if ($msg) {
        print "\n\n";
        print "#################################################################\n";
        print "$msg\n";
        print "#################################################################\n";
    }
    
    print "CMD: $cmd\n";
    if ($PAUSE_EACH_STEP) {
        print STDERR "\n\n-WAITING, PRESS RETURN TO CONTINUE ...";
        my $wait = <STDIN>;
        print STDERR "executing cmd.\n\n";
        
    }
    

    my $time_start = time();
    
    my $ret = system($cmd);
    my $time_end = time();

    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    my $number_minutes = sprintf("%.1f", ($time_end - $time_start) / 60);
    
    print "TIME: $number_minutes min. for $cmd\n";
    

    return;
}


####
sub parse_sample_descriptions {
    my ($read_samples_descr_file, $PARAMS_href) = @_;

    my %samples_descr;
    
    
    open (my $fh, $read_samples_descr_file) or die "Error, cannot open file: $read_samples_descr_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        s/^\s+|\s+$//g;
        chomp;
        my @x = split(/\t/);
        if (/=/ && scalar(@x) == 1) {
            my ($key, $val) = split(/=/, $x[0]);
            $PARAMS_href->{$key} = $val;
        }
        else {
        
            
            my ($condition, $replicate, $reads_left, $reads_right) = @x;
            
            ## remove gzip extension, will convert to gunzipped version later 
            $reads_left =~ s/\.gz$//;
            $reads_right =~ s/\.gz$// if $reads_right;
        
            $samples_descr{$condition}->{$replicate} = [$reads_left, $reads_right];
        }
    }
    
    close $fh;
    
    return(%samples_descr);
    
    
}
