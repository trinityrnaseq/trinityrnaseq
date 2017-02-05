#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::RealBin or die "error, cannot cd to $FindBin::RealBin";



my @files_to_keep = qw (
cleanme.pl 
Makefile
samples.txt

);


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (-f $file && ! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}

`rm -rf ./*__RSEM`;
`rm -rf ./*__kallisto`;
`rm -rf ./*__express`;
`rm -rf ./hsrep1/`;
`rm -rf ./platRep1/`;
`rm -rf ./logRep1/`;
`rm -rf ./dsRep1/`;
`rm -rf ./Trinity.fasta.salmon_quasi.idx/`;

exit(0);
