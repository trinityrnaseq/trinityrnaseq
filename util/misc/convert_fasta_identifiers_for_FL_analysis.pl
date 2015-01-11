#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;


my $usage = "usage: $0 file.fasta\n\n";

my $fasta_file = $ARGV[0] or die $usage;

main: {

    my $fasta_reader = new Fasta_reader($fasta_file);

    my %seen;

    while (my $seq_obj = $fasta_reader->next()) {
        
        my $header = $seq_obj->get_header();
        
        my $seq = $seq_obj->get_sequence();
       
        my ($trans_id, $gene_id, @rest) = split(/\s+/, $header);

                
        if ($seen{$seq}) {
            next;
        }

        $seen{$seq} = 1;

        
        print ">$trans_id;$gene_id\n$seq\n";
    }
    

    exit(0);
}


            
