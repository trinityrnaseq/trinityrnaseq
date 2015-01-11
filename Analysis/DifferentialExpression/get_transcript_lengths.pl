#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;


my $usage = "usage: $0 transcripts.fasta\n\n";

my $trans_fasta = $ARGV[0] or die $usage;

main: {
    
    my $fasta_reader = new Fasta_reader($trans_fasta);
    
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $seq = $seq_obj->get_sequence();

        print join("\t", $acc, length($seq)) . "\n";
    }


    exit(0);
}



    
    
