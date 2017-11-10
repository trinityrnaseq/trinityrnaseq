#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Pipeliner;
use File::Basename;

my $CPU = 2;

my $usage = <<__EOUSAGE;

############################################################
#
# Required:
#
#  --genome <string>         target genome.fasta file
#
#  --samples_file <string>   Trinity samples file
#
# Optional:
#
#  --gtf <string>            annotation in gtf format
#
#  --CPU <int>               multithreading (default: $CPU)
#
#  --nameSorted              sorts bam by read name
#
###########################################################


__EOUSAGE

    ;




my $help_flag;
my $genome_fa;
my $annotation_gtf;
my $samples_file;
my $nameSorted;

&GetOptions ( 'h' => \$help_flag,
              'genome=s' => \$genome_fa,
              'gtf=s' => \$annotation_gtf,
              'samples_file=s' => \$samples_file,
              'CPU=i' => \$CPU,
              'nameSorted' => \$nameSorted,
    );

if ($help_flag) { die $usage; }

unless ($genome_fa && $samples_file) {
    die $usage;
}

$genome_fa = &Pipeliner::ensure_full_path($genome_fa);
$samples_file = &Pipeliner::ensure_full_path($samples_file);
$annotation_gtf = &Pipeliner::ensure_full_path($annotation_gtf) if $annotation_gtf;



main: {
    
    my @read_sets = &parse_samples_file($samples_file);
    
    ## align reads to the mini-genome using hisat2

    ###########################
    # first, build genome index

    my $pipeliner = new Pipeliner(-verbose => 1);
    
    if ($annotation_gtf) {
                
        $pipeliner->add_commands(new Command("hisat2_extract_splice_sites.py $annotation_gtf  > $annotation_gtf.ss",
                                             "$annotation_gtf.ss.ok"));
        
        $pipeliner->add_commands(new Command("hisat2_extract_exons.py $annotation_gtf > $annotation_gtf.exons",
                                             "$annotation_gtf.exons.ok"));
        
        $pipeliner->add_commands(new Command("hisat2-build --exon $annotation_gtf.exons --ss $annotation_gtf.ss -p $CPU $genome_fa $genome_fa",
                                             "$genome_fa.hisat2.build.ok"));
        
        $pipeliner->run();
    }
    else {
                
        $pipeliner->add_commands(new Command("hisat2-build  -p $CPU $genome_fa $genome_fa",
                                             "$genome_fa.hisat2.nogtf.build.ok"));
        
        $pipeliner->run();
    }
    
    #####################
    ## now run alignments
    
    my $aln_checkpoints_dir = "hisat2_aln_chkpts." . basename($genome_fa);
    unless (-d $aln_checkpoints_dir) {
        mkdir($aln_checkpoints_dir) or die "Error, cannot mkdir $aln_checkpoints_dir";
    }

    my $sorted_opt = "";
    my $sorted_token = "c";
    
    if ($nameSorted) {
        $sorted_opt = "-n";
        $sorted_token = "n";
    }
    
    foreach my $read_set_aref (@read_sets) {
        my ($sample_id, $left_fq, $right_fq) = @$read_set_aref;


        my $bamfile = "$sample_id.${sorted_token}Sorted.hisat2." . basename($genome_fa) . ".bam";
        
        $pipeliner->add_commands(new Command("bash -c \"set -eof pipefail; hisat2 --dta -x $genome_fa -p $CPU -1 $left_fq -2 $right_fq | samtools view -Sb -F 4 | samtools sort $sorted_opt -o $bamfile \" ",
                                             "$aln_checkpoints_dir/$bamfile.ok"));
        
    }
    
    $pipeliner->run();
    
    exit(0);
    
}



####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples;
    
    open(my $fh, $samples_file) or die "Error, cannot open file $samples_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($cond, $rep, $fq_a, $fq_b) = @x;

        $fq_a = &Pipeliner::ensure_full_path($fq_a);
        $fq_b = &Pipeliner::ensure_full_path($fq_b) if $fq_b;
        
        push (@samples, [$rep, $fq_a, $fq_b]);
    }
    close $fh;

    return (@samples);
}

