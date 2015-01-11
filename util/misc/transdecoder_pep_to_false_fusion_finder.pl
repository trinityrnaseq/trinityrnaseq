#!/usr/bin/env perl

use strict;
use warnings;

## Just identifies and counts transcripts that are predicted to encode multiple ORFs

# of those that encode an ORF of min length, how many encode multiple ORFs of min length?

my $usage = "\n\nusage: $0 transdecoder.pep [min_prot_length=300]\n\n";

my $transdecoder_pep_file = $ARGV[0] or die $usage;
my $min_prot_length = $ARGV[1] || 300;



main: {

    # example transdecoder pep header
    # m.83 g.83  ORF g.83 m.83 type:complete len:773 (+) CUFF.100.1:5314-7632(+)
    
    my %acc_orf_counter;

    open (my $fh, $transdecoder_pep_file) or die $!;
    while (<$fh>) {
        chomp;
        if (/^>/) {
            /\s(\S+):\d+-\d+\([+-]\)/ or die "Error, cannot extract transcript accession from $_";
            my $trans_acc = $1;

            /\slen:(\d+)\s/ or die "Error, cannot extract len value from $_";
            my $len = $1;

            if ($len >= $min_prot_length) {
                $acc_orf_counter{$trans_acc}++;
            }
        }
    }
    close $fh;

    ## determine number of putative chimeras
    my $num_trans = 0;
    my $num_chims = 0;
    
    foreach my $acc (keys %acc_orf_counter) {

        $num_trans++;

        if ($acc_orf_counter{$acc} > 1) {
            
            $num_chims++;
        
        }

    }

    my $pct_chims = sprintf("%.2f", $num_chims/$num_trans*100);

    print "#trans\tchims\t\%chims\n";
    print join("\t", $num_trans, $num_chims, $pct_chims) . "\n";
    

    exit(0);



}
