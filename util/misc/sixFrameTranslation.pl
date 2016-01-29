#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;
use Nuc_translator;


my $usage = "\nn\\tusage: $0 nucFasta  \n\n";

my $fasta_file = $ARGV[0] or die $usage;

main: {
	my $fasta_reader = new Fasta_reader($fasta_file);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();
		my $header = $seq_obj->get_header();
		
		for my $frame (1..6) {
            my $translation = &translate_sequence($sequence, $frame);
            
            $translation =~ s/(\S{60})/$1\n/g;
            chomp $translation;
            
            print ">$accession.F${frame}_trans $header\n$translation\n";
        }
    }
	
	exit(0);
}



