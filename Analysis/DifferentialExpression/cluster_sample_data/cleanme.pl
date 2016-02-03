#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::RealBin or die "error, cannot cd to $FindBin::RealBin";



my @files_to_keep = qw (cleanme.pl 
                        runMe.sh
                        MLF_ESC_NPC.cuff.genes.fpkm.matrix.gz
                        samples.txt
                        orig.samples.txt
                        );


my %keep = map { + $_ => 1 } @files_to_keep;



foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
