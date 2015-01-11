#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 partitions.gff join_within_dist\n";

my $partitions_gff = $ARGV[0] or die $usage;
my $join_within_dist = $ARGV[1] or die $usage;

my @prev_line;

open (my $fh, $partitions_gff) or die "Error, cannot open file $partitions_gff";
while (<$fh>) {
	chomp;
	
	my @x = split(/\t/);
	
	if (!@prev_line) {
		@prev_line = @x;
		next;
	}
	
	
	if ($x[0] ne $prev_line[0] 
		|| 
		$x[3] - $prev_line[4] > $join_within_dist) {

		print join("\t", @prev_line) . "\n";
		@prev_line = @x;
	}
	else {
		# same contig and within distance
		$prev_line[4] = $x[4]; # update right end.
	}
}

print join("\t", @prev_line) . "\n" if @prev_line;


exit(0);

