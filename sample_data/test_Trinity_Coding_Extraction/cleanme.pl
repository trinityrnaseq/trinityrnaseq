#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 
                        runMe.sh
                        Trinity.fasta.gz
                        __test_using_both_strands.sh
                        __runMe_incl_PfamSearch.sh

                        Makefile
);


my %keep = map { + $_ => 1 } @files_to_keep;



foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


`rm -rf ./transdecoder.tmp.*`;

exit(0);
