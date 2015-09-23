#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $usage = "usage: $0 (RSEM|eXpress|kallisto)\n\n";

my $method = $ARGV[0] or die $usage;
unless ($method =~ /^(RSEM|eXpress|kallisto)$/) {
    die $usage;
}

my $utildir = "$FindBin::Bin/../../util";

my $samples_file = "samples_n_reads_decribed.txt";
my $trinity_fasta = "../test_DATA/Trinity.fasta";

main: {
    
    my @samples;
    {
        open (my $fh, $samples_file) or die $!;
        while (<$fh>) {
            chomp;
            my ($sample_name, $left_fq, $right_fq) = @_;
            push (@samples, [$sample_name, $left_fq, $right_fq]);
        }
        close $fh;
    }

    my @trans_matrices;
    my @gene_matrices;
    
    foreach my $sample (@samples) {
        my ($sample_name, $left_fq, $right_fq) = @$sample;

        my $cmd = "$utildir/align_and_estimate_abundance.pl --transcripts $trinity_fasta "
            . " --left $left_fq --right $right_fq --seqType fq ";
        
        if ($method eq 'RSEM') {
            $cmd .= " --est_method RSEM --output_dir RSEM-$sample_name --aln_method bowtie ";

        }
        elsif ($method eq 'eXpress') {
            $cmd .= " --est_method eXpress --output_dir eXpress-$sample_name --aln_method bowtie2 ";
            

        }
        elsif ($method eq 'kallisto') {
            $cmd .= " --est_method kallisto --output_dir kallisto-$sample_name ";


        }
        else {
            # shouldn't ever get here.
            die "error - method $method not recognized";
        }
    }

    ## generate matrices.
    

    exit(0);
}

