#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;
use Nuc_translator;

my $usage = "\n\nusage: $0 transcripts.fa\n\n";

my $transcripts_file = $ARGV[0] or die $usage;


main: {

    my $fasta_reader = new Fasta_reader($transcripts_file);

    while (my $seq_obj = $fasta_reader->next()) {
        
        my $acc = $seq_obj->get_accession();
        
        my $seq = $seq_obj->get_sequence();

        print ">$acc\n$seq\n";

        my $revc_seq = &reverse_complement($seq);

        print ">$acc-ANTI\n$revc_seq\n";
        
    }
    

    exit(0);

}
