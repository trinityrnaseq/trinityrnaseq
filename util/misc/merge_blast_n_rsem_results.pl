#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 rsem.out blast.outfmt6 [transcripts.fasta]\n\n";

my $rsem_out = $ARGV[0] or die $usage;
my $blast_out = $ARGV[1] or die $usage;
my $transcripts_fasta = $ARGV[2];

main: {

    my %rsem_text;
    my $rsem_header;
    {
        open (my $fh, $rsem_out) or die $!;
        $rsem_header = <$fh>;
        chomp $rsem_header;
        while (<$fh>) {
            chomp;
            my $line = $_;
            my @x = split(/\t/);
            my ($gene, $trans_list) = ($x[0], $x[1]);
            foreach my $ele ($gene, split(/,/, $trans_list)) {
                $rsem_text{$ele} = $line;
            }
        }
        close $fh;
    }

    my %trans_seqs;
    if ($transcripts_fasta) {
        my $fasta_reader = new Fasta_reader($transcripts_fasta);
        %trans_seqs = $fasta_reader->retrieve_all_seqs_hash($transcripts_fasta);
    }
    

    open (my $fh, $blast_out) or die $!;
    my $header = <$fh>;
    print "$rsem_header\t$header";
    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $acc = $x[0];
        my $rsem_line = $rsem_text{$acc} or die "Error, no rsem text for $acc";
        print "$rsem_line\t$line";
        if (my $seq = $trans_seqs{$acc}) {
            print "\t$seq";
        }
        print "\n";
    }
    close $fh;
    
    exit(0);
}
