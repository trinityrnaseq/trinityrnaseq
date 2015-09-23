#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Basename;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;


my $usage = "usage: $0 (RSEM|eXpress|kallisto)\n\n";

my $method = $ARGV[0] or die $usage;
unless ($method =~ /^(RSEM|eXpress|kallisto)$/) {
    die $usage;
}

my $utildir = "$FindBin::Bin/../../util";

my $samples_file = "samples.txt";
my $trinity_fasta = "../test_DATA/Trinity.fasta";

main: {

    &process_cmd("ln -sf $trinity_fasta .");

    $trinity_fasta = basename($trinity_fasta);
    
    my @samples;
    {
        open (my $fh, $samples_file) or die $!;
        while (<$fh>) {
            chomp;
            my ($sample_name, $left_fq, $right_fq) = split(/\s+/);
            $left_fq = &ensure_full_path($left_fq);
            $right_fq = &ensure_full_path($right_fq);
            push (@samples, [$sample_name, $left_fq, $right_fq]);
        }
        close $fh;
    }

    my @trans_results;
    my @gene_results;
    
    foreach my $sample (@samples) {
        my ($sample_name, $left_fq, $right_fq) = @$sample;

        my $cmd = "$utildir/align_and_estimate_abundance.pl --transcripts $trinity_fasta --prep_reference "
            . " --left $left_fq --right $right_fq --seqType fq --trinity_mode ";
        
        my $outdir = "$method-$sample_name";

        if ($method eq 'RSEM') {
            $cmd .= " --est_method RSEM --output_dir $outdir --aln_method bowtie ";
            
            push (@trans_results, "$outdir/RSEM.isoforms.results");
            push (@gene_results, "$outdir/RSEM.genes.results");
            
        }
        elsif ($method eq 'eXpress') {
            $cmd .= " --est_method eXpress --output_dir $outdir --aln_method bowtie2 ";
            push (@trans_results, "$outdir/results.xprs");
            push (@gene_results, "$outdir/results.xprs.genes");
        }
        elsif ($method eq 'kallisto') {
            $cmd .= " --est_method kallisto --output_dir kallisto-$sample_name ";
            push (@trans_results, "$outdir/abundance.tsv");
            push (@gene_results, "$outdir/abundance.tsv.genes");
            
        }
        else {
            # shouldn't ever get here.
            die "error - method $method not recognized";
        }

        &process_cmd($cmd);


    }

    ## generate matrices.
    my $cmd = "$utildir/abundance_estimates_to_matrix.pl --est_method $method --out_prefix $method-trans --name_sample_by_basedir @trans_results";
    &process_cmd($cmd);
    
    $cmd = "$utildir/abundance_estimates_to_matrix.pl --est_method $method --out_prefix $method-gene --name_sample_by_basedir @gene_results";
    &process_cmd($cmd);
    
    

    exit(0);
}

