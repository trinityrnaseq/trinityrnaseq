#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib "$FindBin::Bin/../../../PerlLib";
use Pipeliner;
use File::Basename;

my $CPU = 2;

my $top_genes_plot = 50;

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
#  --SS_lib_type <string>            strand-specific library type 'RF|FR|R|F'
#
#  --CPU <int>                       default: $CPU
#
#  --top_genes_plot <int>            default: $top_genes_plot
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
my $SS_lib_type = "";

&GetOptions ( 'h' => \$help_flag,
              'genes_fasta=s' => \$genes_fasta_file,
              'genes_gtf=s' => \$genes_gtf_file,
              'samples_file=s' => \$samples_file,
              'out_prefix=s' => \$out_prefix,
              'CPU=i' => \$CPU,
              'aligner=s' => \$aligner,
              'SS_lib_type=s' => \$SS_lib_type,
              'top_genes_plot=i' => \$top_genes_plot,
              
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

    
    my $chkpts_dir = "dexseq_chkpts";
    unless(-d $chkpts_dir) {
        mkdir $chkpts_dir or die "Error, cannot mkdir $chkpts_dir";
    }
    
    my $analysis_token = "$chkpts_dir/" . basename($genes_gtf_file) . ".$aligner";
    
    ## flatten the gtf file
    my $cmd = "$TRINITY_HOME/trinity-plugins/DEXseq_util/dexseq_prepare_annotation.py $genes_gtf_file $genes_gtf_file.dexseq.gff";
    $pipeliner->add_commands(new Command($cmd, "$chkpts_dir/" . basename($genes_gtf_file) . ".flatten_gtf.ok"));

    if ($aligner =~ /HISAT2/i) {
        $cmd = "$TRINITY_HOME/util/misc/run_HISAT2_via_samples_file.pl --genome $genes_fasta_file --gtf $genes_gtf_file --samples_file $samples_file --CPU $CPU --nameSorted ";
        $pipeliner->add_commands(new Command($cmd, "$analysis_token.align.ok"));
        
        $pipeliner->run();
    }
    elsif ($aligner =~ /STAR/i) {
        $cmd = "$TRINITY_HOME/util/misc/run_STAR_via_samples_file.pl --genome $genes_fasta_file --gtf $genes_gtf_file --samples_file $samples_file --CPU $CPU --nameSorted ";
        $pipeliner->add_commands(new Command($cmd, "$analysis_token.align.ok"));
        $pipeliner->run();
    }
    else {
        die "Error, dont recognize aligner: $aligner";
    }

    my $genes_file_basename = basename($genes_fasta_file);
    
    my @counts_files;
    ## process each of the replicates
    foreach my $sample_info_aref (@samples_info) {
        my ($condition, $replicate, $left_fq, $right_fq) = @$sample_info_aref;
        
        my $bam_file = "$replicate.nSorted.__METHOD__.__GENES__.bam";
        $bam_file =~ s/__METHOD__/$aligner/;
        $bam_file =~ s/__GENES__/$genes_file_basename/;
        
        unless (-s $bam_file) {
            die "Error, cannot locate bam file: $bam_file";
        }
        
        # quant
        my $cmd = "featureCounts -T $CPU ";



        # paired or unpaired:
        if ($right_fq) {
            $cmd .= " -p -B ";
        }

        ## strand-specific options:
        if ($SS_lib_type) {
            if ($SS_lib_type =~ /^F/) {
                $cmd .= " -s 1";
            }
            elsif ($SS_lib_type =~ /^R/) {
                $cmd .= " -s 2";
            }
            else {
                die "Error, not recognizing strand-specific lib type [$SS_lib_type]";
            }
        }
        else {
            # unstranded by default
        }
        
        $cmd .= " -O --fraction -M -t exonic_part -g exonic_gene_part_number -a $genes_gtf_file.dexseq.gff -f -o $bam_file.fc  $bam_file";
        
        $pipeliner->add_commands(new Command($cmd, "$analysis_token.$bam_file.fc.ok"));

        
        $pipeliner->add_commands(new Command("$FindBin::Bin/util/reformat_featureCounts.pl $bam_file.fc > $bam_file.counts",
                                             "$analysis_token.$bam_file.counts.ok"));
        
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
            my ($condition, $replicate_name, $left_fq, $right_fq) = @$sample_info_aref;
            my $counts_file = $counts_files[$i];

            print $ofh join("\t", $replicate_name, $condition, $counts_file) . "\n";
        }

        close $ofh;
    }
    
    my $dexseq_rscript = "$out_prefix.dexseq.Rscript";
    if ($out_prefix !~ /dexseq/i) {
        $out_prefix .= ".dexseq";
    }
    
    {
        
        open (my $ofh, ">$dexseq_rscript") or die "Error, cannot write to $dexseq_rscript";
        print $ofh "library(DEXSeq)\n";
        print $ofh "samples_info = read.table(\"$samples_table_file\", header=T, row.names=1)\n";
        print $ofh "dxd = DEXSeqDataSetFromHTSeq(as.vector(samples_info\$counts_filename), sampleData=samples_info, design = ~ sample + exon + condition:exon, flattenedfile=\"$genes_gtf_file.dexseq.gff\")\n";
        print $ofh "dxd = estimateSizeFactors( dxd )\n";
        print $ofh "dxd = estimateDispersions( dxd )\n";
        print $ofh "plotDispEsts( dxd )\n";
        print $ofh "dxd = testForDEU( dxd )\n";
        print $ofh "dxd = estimateExonFoldChanges( dxd, fitExpToVar=\"condition\")\n";
        print $ofh "dxr1 = DEXSeqResults( dxd )\n";
        print $ofh "dxr1.sorted = dxr1[order(dxr1\$padj),]\n";
        print $ofh "save(list = ls(all=TRUE), file = \"$out_prefix.Rdata\")\n";
        print $ofh "write.table(dxr1.sorted, file=\"$out_prefix.results.dat\", quote=F, sep=\"\t\")\n";
        
        print $ofh "pdf(\"$out_prefix.pdf\")\n";
        
        print $ofh "top_genes = unique(dxr1.sorted\$groupID[dxr1.sorted\$padj < 0.1 & ! is.na(dxr1.sorted\$padj)])\n";
        print $ofh "top_genes = top_genes[1:min($top_genes_plot, length(top_genes))]\n";
        print $ofh "message(\"Top $top_genes_plot genes: (\", paste(top_genes, collapse=','), \")\")\n"; 
        print $ofh "for (gene in top_genes) { \n"
            . "    plotDEXSeq( dxr1 , gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , expression=FALSE, norCounts=TRUE, splicing=TRUE, displayTranscripts=TRUE)\n"
            . "}\n";
                        
        #print $ofh "plotMA( dxr1, cex=0.8 )\n";
        close $ofh;
        
    }
    
    $cmd = "R --no-save --no-restore --no-site-file --no-init-file -q < $dexseq_rscript";
    $pipeliner->add_commands(new Command($cmd, "$analysis_token.dexseq.$top_genes_plot.ok"));
    
    $pipeliner->run();

    
}

####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples_info;

    my %seen;
    
    open (my $fh, $samples_file) or die "Error, cannot open file: $samples_file";
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        chomp;
        my $line = $_;
        my @x = split(/\s+/);
        my $condition = $x[0];
        my $replicate = $x[1];
        my $left_fq = $x[2];
        my $right_fq = $x[3]; # only if paired-end

        unless ($left_fq) {
            die "Error, cannot find required fields in sample file line: $line";
        }
        unless (-s $left_fq) {
            die "Error, cannot locate left_fq: $left_fq as specified in samples file line: $line ";
        }
        if ($right_fq && ! -s $right_fq) {
            die "Error, cannot locate right_fq: $right_fq as specified in samples file line: $line";
        }
        
        if ($seen{$replicate}) {
            die "Error, replicate name: [$replicate] must be unique among all replicate names.  Please update your samples file";
        }
        $seen{$replicate} = 1;
        
        push (@samples_info, [$condition, $replicate, $left_fq, $right_fq]);
    }
    close $fh;

    return(@samples_info);
}
