#!/usr/bin/env perl

use strict;

my $delimeter = ($ARGV[0] eq '-s') ? '\s+' : '\t';

while (<STDIN>) {
    
	unless (/\w/) {
		print "0\n\n";
		next;
	}

	chomp;
    my $x = 0;
    my @x = split (/$delimeter/);
    foreach my $word (@x) {
		print "$x\t$word\n";
		$x++;
    }
    print "\n";
    #last;
}
