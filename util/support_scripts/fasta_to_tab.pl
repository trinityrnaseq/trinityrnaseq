#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 [multiFastaFile] [NO_FULL_HEADER_FLAG=0]\n\n";

my $input = $ARGV[0] || *STDIN{IO};
my $NO_FULL_HEADER_FLAG = $ARGV[1] || 0;

unless (-f $input || ref $input eq 'IO::Handle') {
	die "Error, input not established.";
}

main: {
	
	my $fasta_reader = new Fasta_reader($input);
	
	while (my $seq_obj = $fasta_reader->next()) {
		my $sequence = $seq_obj->get_sequence();
		my $header = $seq_obj->get_header();

        if ($NO_FULL_HEADER_FLAG) {
            $header = $seq_obj->get_accession();
        }
        
		print "$header\t$sequence\n";
	}
    
	exit(0);
}



