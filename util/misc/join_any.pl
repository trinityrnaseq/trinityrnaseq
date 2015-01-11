#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage: $0 token_list_file  file_to_join [-v]\n\n";

my $token_list_file = $ARGV[0] or die $usage;
my $file_to_join = $ARGV[1] or die $usage;
my $invert_selection = $ARGV[2] || 0;

my %tokens;
{ 
	open (my $fh, $token_list_file) or die "Error, cannot open file $token_list_file";
	while (<$fh>) {
		while (/(\S+)/g) {
			$tokens{$1} = 1;
		}
	}
	close $fh;
}


open (my $fh, $file_to_join) or die "Error, cannot open file $file_to_join ";
while (<$fh>) {
	my $line = $_;
	chomp;
	my @x = split (/\s+/);
	my $found_token = 0;
	foreach my $ele (@x) {
		if ($tokens{$ele}) {
			$found_token = 1;
			last;
		}
	}

	if ($found_token && !$invert_selection) {
		print $line;
	}
	elsif ($invert_selection && !$found_token) {
		print $line;
	}
	
}
close $fh;


exit(0);


		
