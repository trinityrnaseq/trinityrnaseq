#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "\n\n\tusage: $0 < list of audit.txt files from stdin \n\n\n";

my $MIN_SEQ_LEN = 1000;



my $total_genes = 0;
my $total_genes_reco = 0;

my $total_refseq_trans = 0;
my $total_iso_reco_count = 0;
my $total_extra = 0;

print join("\t", "gene_id", "reco_gene_flag", "num_refseqs", "num_FL_trin", "num_extra_trin") . "\n";

my $counter = 0;
while (<>) {
    chomp;
    my $dir = dirname($_);

    $counter++;
    
    my $gene_id = basename($dir);

    my $trin_fasta_file = "$dir/trinity_out_dir.Trinity.fasta";

    if (! -e $trin_fasta_file) {
        print STDERR "ERROR: Cannot locate file: $trin_fasta_file\n";
        next;
    }
    
    my $fasta_reader = new Fasta_reader($trin_fasta_file);
    my %trin_seqs = $fasta_reader->retrieve_all_seqs_hash();
    
    my %trin_lens = &get_seq_lens(%trin_seqs);

    my $FL_reco_file = "$dir/FL.test.pslx.maps";
    my %reco_refseq;
    my %reco_trin_to_refseq = &parse_reco($FL_reco_file, \%reco_refseq);
    my $num_FL_trin = scalar(keys %reco_trin_to_refseq);

    my $refseqs_fa = "$dir/refseqs.fa";
    my @refseq_accs = &get_accs($refseqs_fa);

    my $num_refseqs = scalar(@refseq_accs);
    
    my @failed_reco_refseqs;

    foreach my $acc (@refseq_accs) {
        unless ($reco_refseq{$acc}) {
            push (@failed_reco_refseqs, $acc);
        }
    }
    my $num_failed_reco_refseqs = scalar(@failed_reco_refseqs);
    
    my @extra_trin_accs;
    foreach my $trin_acc (keys %trin_seqs) {
        if (! exists $reco_trin_to_refseq{$trin_acc}) {
            if ($trin_lens{$trin_acc} >= $MIN_SEQ_LEN) {
                push (@extra_trin_accs, $trin_acc);
            }
        }
    }

    my $num_extra_trin = scalar(@extra_trin_accs);

    my $reco_gene_flag = ($num_failed_reco_refseqs == 0) ? "YES" : "NO";
    
    print join("\t", $gene_id, $reco_gene_flag, $num_refseqs, $num_FL_trin, $num_extra_trin) . "\n";
    
    $total_genes++;
    if ($reco_gene_flag eq "YES") {
        $total_genes_reco++;
    }
    $total_refseq_trans += $num_refseqs;
    $total_iso_reco_count += $num_FL_trin;
    $total_extra += $num_extra_trin;

}

if ($counter > 1) {
    print "\n\n";
    print join("\t", "Total_Genes", "Total_Genes_Reco", "Total_RefTrans", "Total_RefTransReco", "Total_extra_trans") . "\n";
    print join("\t", $total_genes, $total_genes_reco, $total_refseq_trans, $total_iso_reco_count, $total_extra) . "\n";
}


exit(0);


####
sub get_accs {
    my ($fasta_file) = @_;

    my @accs;

    open (my $fh, $fasta_file) or die $!;
    while (<$fh>) {
        if (/^>(\S+)/) {
            push (@accs, $1);
        }
    }

    return(@accs);
}


####
sub parse_reco {
    my ($reco_file, $reco_refseq_href) = @_;

    my %trin_to_reco_acc;
    
    open (my $fh, $reco_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($trans_acc, $trinity_contigs) = split(/\t/);
        
        my @trin_contigs = split(/,/, $trinity_contigs);
        foreach my $trin (@trin_contigs) {
            
            $trin_to_reco_acc{$trin}->{$trans_acc} = 1;
            
            $reco_refseq_href->{$trans_acc} = 1;
        }
    }
    close $fh;
    
    return(%trin_to_reco_acc);
    
}


####
sub get_seq_lens {
    my (%trin_seqs) = @_;

    my %lens;

    foreach my $acc (keys %trin_seqs) {
        my $seq = $trin_seqs{$acc};
        my $seqlen = length($seq);
        
        $lens{$acc} = $seqlen;
    }

    return(%lens);

}
    
