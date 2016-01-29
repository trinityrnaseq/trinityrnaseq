#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use FindBin;
use Cwd;

######################################################
## Set to base directory of the Trinity installation:
my $BASEDIR = "$FindBin::Bin/../";
######################################################

my $usage = <<__EOUSAGE__;


##########################################################################################################
#
#  Required:
#
#  --trinity_fasta <string>        Trinity fasta file
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
my $trinity_fasta_file;

&GetOptions ( 'h' => \$help_flag,
              'trinity_fasta=s' => \$trinity_fasta_file,
              'samples_file=s' => \$read_samples_descr_file,
              'I' => \$PAUSE_EACH_STEP,
              
);


if ($help_flag) {
    die $usage;
}

unless ($trinity_fasta_file && $read_samples_descr_file) {
    die $usage;
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

#############################################################
## Run RSEM, compute abundance estimates for each replicate.
#############################################################



my @trans_rsem_files;
my @genes_rsem_files;

foreach my $condition (@conditions) {

    my $replicates_href = $conditions_to_read_info{$condition};

    my @replicates = keys %$replicates_href;
    foreach my $replicate (@replicates) {
        my ($left_fq_file, $right_fq_file) = @{$replicates_href->{$replicate}};
        
        ## run RSEM
        my $seqType = $PARAMS{"--seqType"} or die "Error, --seqType not specified";
        my $cmd = "$BASEDIR/util/align_and_estimate_abundance.pl "
            . " --transcripts $trinity_fasta_file "
            . " --seqType $seqType "
            . " --prep_reference "
            . " --output_prefix $replicate "
            . " --aln_method bowtie --est_method RSEM "
            . " --trinity_mode "
            
            ;
        
        if ($left_fq_file && $right_fq_file) {
            $cmd .= " --left $left_fq_file --right $right_fq_file ";
        }
        elsif ($left_fq_file) {
            $cmd .= " --single $left_fq_file ";
        }
        else {
            die "Error, no left or right read files";
        }
        if (my $SS_lib_type = $PARAMS{'--SS_lib_type'}) {
            $cmd .= " --SS_lib_type $SS_lib_type ";
        }
        if (my $cpu = $PARAMS{'--CPU'}) {
            $cmd .= " --thread_count $cpu ";
        }
        
        &process_cmd($cmd,
                     "Running RSEM to perform abundance estimation on sample $replicate"
                     ) unless (-s "$replicate.isoforms.results");
        
            
        push (@trans_rsem_files, "$replicate.isoforms.results");
        push (@genes_rsem_files, "$replicate.genes.results");
    }
    
}


######################################
# Make count matrices for DE analysis
######################################

# make iso matrix
my $cmd = "$BASEDIR/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix Trinity_trans " . join (" " , @trans_rsem_files);

&process_cmd($cmd,
             "Generating transcript (isoform) count matrix file."
             ) unless (-s "Trinity_trans.counts.matrix");


# make the 'gene' matrix
$cmd = "$BASEDIR/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix Trinity_genes " . join (" " , @genes_rsem_files);

&process_cmd($cmd,
             "Generating Trinity component (gene) count matrix file"
    ) unless (-s "Trinity_genes.counts.matrix");




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
    
    
    open (my $fh, $read_samples_descr_file) or die $!;
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
