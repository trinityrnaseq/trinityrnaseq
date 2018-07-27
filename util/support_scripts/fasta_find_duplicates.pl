#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;


my $usage = "usage: $0 seqs.fasta\n\n";

my $file = $ARGV[0] or die $usage;

my %seq_to_header;

my $fasta_reader = new Fasta_reader($file);
while (my $seq_obj = $fasta_reader->next()) {
    
    my $sequence = $seq_obj->get_sequence();
    my $header = $seq_obj->get_header();
    
    if (exists $seq_to_header{$sequence}) {
        push (@{$seq_to_header{$sequence}}, $header);
    }
    else {
        $seq_to_header{$sequence} = [$header];
    }
}

my $found_dups_flag = 0;
foreach my $sequence (keys %seq_to_header) {
    my @dups  = @{$seq_to_header{$sequence}};

    if (scalar(@dups) > 1) {
        print "# Repeated seqs found:\n";
        print join("\n", @dups) . "\n";
        print "$sequence\n\n";
        $found_dups_flag = 1;
    }
    
}


exit($found_dups_flag);


