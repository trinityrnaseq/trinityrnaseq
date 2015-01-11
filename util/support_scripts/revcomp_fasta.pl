#!/usr/bin/env perl

use strict;
use warnings;

local $/ = ">";
while(<>){
	chomp;
	my $lineSepPos = index($_, "\n");
	my $header = substr($_,0,$lineSepPos);
	if($header){
		print(">", $header, "\n");
		my $sequence = reverse(substr($_,$lineSepPos));
		# see http://shootout.alioth.debian.org/u32/performance.php?test=revcomp#about
		# for ambiguity codes and translation
		$sequence =~ tr/ACGTUMRWSYKVHDBNacgtumrwsykvhdbn\n/TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn/d;
		for(my $pos = 0; $pos < length($sequence);$pos += 60){
			print(substr($sequence, $pos, 60),"\n");
		}
	}
}
