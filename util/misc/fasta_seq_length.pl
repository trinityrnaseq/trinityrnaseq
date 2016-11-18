#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 fastaFile\n\n";

my $file = $ARGV[0] or die $usage;

my $fasta_reader = new Fasta_reader($file);

print join("\t", "#fasta_entry", "length") . "\n";
while (my $seq_obj = $fasta_reader->next()) {
    my $sequence = $seq_obj->get_sequence();
    my $accession = $seq_obj->get_accession();
    
    print join("\t", $accession, length($sequence)) . "\n";
}

exit(0);

