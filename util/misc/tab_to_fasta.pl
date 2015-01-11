#!/usr/bin/env perl

use strict;
use warnings;

main: {
	while (<STDIN>) {
		unless (/\w/) { next; }
		chomp;
		my @x = split (/\t/);
		my $seq = pop @x;

		# make fasta
		$seq =~ s/(\S{60})/$1\n/g;
		chomp $seq;

		print ">" . join (" ", @x) . "\n$seq\n";
		
	}

	exit(0);
}


