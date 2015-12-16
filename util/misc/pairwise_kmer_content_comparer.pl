#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 file.fasta [kmer_length=25]\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $kmer_size = $ARGV[1] || 25;


main: {

    my $fasta_reader = new Fasta_reader($fasta_file);
    my %trans_seqs = $fasta_reader->retrieve_all_seqs_hash();

    my %acc_to_kmers = &get_kmers(\%trans_seqs);

    my @accs = keys %acc_to_kmers;

    
    for (my $i = 0; $i < $#accs; $i++) {
        
        my $acc_i = $accs[$i];
        my $kmers_href_i = $acc_to_kmers{$acc_i};


        for (my $j = $i + 1; $j <= $#accs; $j++) {
            
            my $acc_j = $accs[$j];
            my $kmers_href_j = $acc_to_kmers{$acc_j};

            if (&have_kmer_in_common($kmers_href_i, $kmers_href_j)) {

                print join("\t", $acc_i, $acc_j) . "\n";
            }
        }

    }

    exit(0);
    
}

####
sub get_kmers {
    my ($trans_seqs_href) = @_;

    my %acc_to_kmers;

    my $counter = 0;
    foreach my $acc (keys %$trans_seqs_href) {
        $counter++;
        print STDERR "-getting kmers for $counter\n";


        my $trans_seq = $trans_seqs_href->{$acc};
        
        my %kmers = &extract_kmers($trans_seq);
        $acc_to_kmers{$acc} = \%kmers;
    }
    
    return(%acc_to_kmers);
}

####
sub extract_kmers {
    my ($seq) = @_;

    my %kmers;

    for (my $i = 0; $i <= length($seq) - $kmer_size; $i++) {
        my $kmer = substr($seq, $i, $kmer_size);
        
        $kmers{$kmer} = 1;
    }

    return(%kmers);
}

####
sub have_kmer_in_common {
    my ($kmers_i, $kmers_j) = @_;

    foreach my $kmer (keys %$kmers_i) {

        if (exists $kmers_j->{$kmer}) {
            return(1);
        }
    }

    return(0);
}

