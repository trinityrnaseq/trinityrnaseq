#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib "$FindBin::Bin/../../../PerlLib";
use Pipeliner;


my $CPU = 2;

my $usage = <<__EOUSAGE__;

#################################################################
#
#  Required:
#
#  --genes_fasta <string>     Trinity genes fasta files
#
#  --genes_gtf <string>       Trinity genes gtf file
#
#  --samples_file <string>            Trinity samples file
#
#  --aligner <string>         aligner to use: STAR|HISAT2
#
#  Optional:
#
#  --out_prefix <string>             default: 'dexseq'
#
#  --CPU <int>                       default: $CPU
#
################################################################


__EOUSAGE__

    ;

#  --aligner <string>                hisat2|gsnap|STAR


my $help_flag;
my $genes_fasta_file;
my $genes_gtf_file;
my $samples_file;
my $out_prefix = "dexseq";
my $aligner;

&GetOptions ( 'h' => \$help_flag,
              'genes_fasta=s' => \$genes_fasta_file,
              'genes_gtf=s' => \$genes_gtf_file,
              'samples_file=s' => \$samples_file,
              'out_prefix=s' => \$out_prefix,
              'CPU=i' => \$CPU,
              'aligner=s' => \$aligner,
    );

if ($help_flag) {
    die $usage;
}

unless ($genes_fasta_file && $genes_gtf_file && $samples_file && $aligner) {
    die $usage;
}
unless ($aligner =~ /^(STAR|HISAT2)$/i) {
    die "Error, dont recognize aligner [$aligner]";
}

$aligner = lc $aligner;

my $TRINITY_HOME = "$FindBin::Bin/../../..";

main: {

    my @samples_info = &parse_samples_file($samples_file);

    my $pipeliner = new Pipeliner(-verbose => 2);

    ## flatten the gtf file
    my $cmd = "$TRINITY_HOME/trinity-plugins/DEXseq_util/dexseq_prepare_annotation.py $genes_gtf_file $genes_gtf_file.dexseq.gff";
    $pipeliner->add_commands(new Command($cmd, "flatten_gtf.ok"));
    

    if ($aligner =~ /HISAT2/i) {
        $cmd = "$TRINITY_HOME/util/misc/run_HISAT2_via_samples_file.pl --genome $genes_fasta_file --gtf $genes_gtf_file --samples_file $samples_file --CPU $CPU ";
        $pipeliner->add_commands(new Command($cmd, "hisat2.ok"));
        
        $pipeliner->run();
    }
    elsif ($aligner =~ /STAR/i) {
        $cmd = "$TRINITY_HOME/util/misc/run_STAR_via_samples_file.pl --genome $genes_fasta_file --gtf $genes_gtf_file --samples_file $samples_file --CPU $CPU ";
        $pipeliner->add_commands(new Command($cmd, "star.ok"));
        $pipeliner->run();
    }
    else {
        die "Error, dont recognize aligner: $aligner";
    }
        
    my @counts_files;
    ## process each of the replicates
    foreach my $sample_info_aref (@samples_info) {
        my ($condition, $replicate) = @$sample_info_aref;
        
        my $bam_file = "$replicate.cSorted.__METHOD__.bam";
        $bam_file =~ s/__METHOD__/$aligner/;
        
        unless (-s $bam_file) {
            die "Error, cannot locate bam file: $bam_file";
        }
        
        # convert to sam
        my $cmd = "samtools view $bam_file > $bam_file.sam";
        $pipeliner->add_commands(new Command($cmd, "$bam_file.sam.ok"));

        # quant
        $cmd = "$TRINITY_HOME/trinity-plugins/DEXseq_util/dexseq_count.py $genes_gtf_file.dexseq.gff $bam_file.sam $bam_file.counts";
        $pipeliner->add_commands(new Command($cmd, "$bam_file.counts.ok"));
        
        push (@counts_files, "$bam_file.counts");
        
    }

    #############
    ## run DEXseq
    
    # first write a samples table
    my $samples_table_file = "$out_prefix.sample_table";
    {
        open (my $ofh, ">$samples_table_file") or die "Error, cannot write to $samples_table_file";
        
        print $ofh "\t" . join("\t", "condition", "counts_filename") . "\n";
        for (my $i = 0; $i <= $#samples_info; $i++) {
            my $sample_info_aref = $samples_info[$i];
            my ($condition, $replicate_name) = @$sample_info_aref;
            my $counts_file = $counts_files[$i];

            print $ofh join("\t", $replicate_name, $condition, $counts_file) . "\n";
        }

        close $ofh;
    }
    
    my $dexseq_rscript = "$out_prefix.Rscript";
    {
        
        open (my $ofh, ">$dexseq_rscript") or die "Error, cannot write to $dexseq_rscript";
        print $ofh "library(DEXSeq)\n";
        print $ofh "samples_info = read.table(\"$samples_table_file\", header=T, row.names=1)\n";
        print $ofh "dxd = DEXSeqDataSetFromHTSeq(as.vector(samples_info\$counts_filename), sampleData=samples_info, design = ~ sample + exon + condition:exon, flattenedfile=\"$genes_gtf_file.dexseq.gff\")\n";
        print $ofh "pdf(\"$out_prefix.pdf\")\n";
        print $ofh "dxd = estimateSizeFactors( dxd )\n";
        print $ofh "dxd = estimateDispersions( dxd )\n";
        print $ofh "plotDispEsts( dxd )\n";
        print $ofh "dxd = testForDEU( dxd )\n";
        print $ofh "dxd = estimateExonFoldChanges( dxd, fitExpToVar=\"condition\")\n";
        print $ofh "dxr1 = DEXSeqResults( dxd )\n";
        print $ofh "write.table(dxr1, file=\"$out_prefix.results.dat\", quote=F, sep=\"\t\")\n";
        print $ofh "plotMA( dxr1, cex=0.8 )\n";
        close $ofh;
        
    }
    
    $cmd = "R --vanilla -q < $dexseq_rscript";
    $pipeliner->add_commands(new Command($cmd, "dexseq.ok"));

    $pipeliner->run();

    
}

####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples_info;
    
    open (my $fh, $samples_file) or die "Error, cannot open file: $samples_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $condition = $x[0];
        my $replicate = $x[1];

        push (@samples_info, [$condition, $replicate]);
    }
    close $fh;

    return(@samples_info);
}
