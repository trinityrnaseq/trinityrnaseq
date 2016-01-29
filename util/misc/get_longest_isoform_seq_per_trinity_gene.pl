#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;

print STDERR "\n\n\tNOTE - longest transcript isn't always the best transcript!... consider filtering based on relative expression support ... \n\n";

my $usage = "usage: $0 Trinity.fasta\n\n";

my $trin_fasta = $ARGV[0] or die $usage;

main: {
    my %gene_to_longest_transcript;
    
    my $fasta_reader = new Fasta_reader($trin_fasta);
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $header = $seq_obj->get_header();
        
        my $gene_id;
        
        if ($acc =~ /^(.*comp\d+_c\d+)_seq/) {
            $gene_id = $1;
        }
        elsif ($acc =~ /^(.*c\d+_g\d+)_i/) {
            $gene_id = $1;
        }

        unless ($gene_id) {
            die "Error, cannot parse gene identifier from acc: $acc";
        }

        my $sequence  = $seq_obj->get_sequence();
        
        if ( (! exists $gene_to_longest_transcript{$gene_id}) 
             ||
             $gene_to_longest_transcript{$gene_id}->{length} < length($sequence)) {

            $gene_to_longest_transcript{$gene_id} = { length => length($sequence), 
                                                      acc => $acc,
                                                      sequence => $sequence,
                                                      header => $header,
            };
        }

        
    }

    foreach my $longest_seq ( reverse sort {$a->{length}<=>$b->{length}} values %gene_to_longest_transcript) {

        print ">" . $longest_seq->{header} . "\n" . $longest_seq->{sequence} . "\n";

    }

    print STDERR "\n\nok.  Done.\n\n";

    exit(0);
}
