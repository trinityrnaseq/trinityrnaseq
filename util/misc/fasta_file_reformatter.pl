#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 fasta\n";

my $fasta_file = $ARGV[0] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);

while (my $seq_obj = $fasta_reader->next()) {
	
	my $header = $seq_obj->get_header();
	my $seq = $seq_obj->get_sequence();

	$seq =~ s/(\S{60})/$1\n/g;
	
    chomp $seq;
	print ">$header\n$seq\n";
	
}

exit(0);

