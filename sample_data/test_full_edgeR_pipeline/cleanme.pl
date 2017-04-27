#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::RealBin or die "error, cannot cd to $FindBin::RealBin";



my @files_to_keep = qw (cleanme.pl 
                        runMe.sh
                        samples_n_reads_decribed.txt
                        
                        Makefile
                        );


my %keep = map { + $_ => 1 } @files_to_keep;

foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}

`rm -rf ./plat_rep1`;
`rm -rf ./hs_rep1`;
`rm -rf ./log_rep1`;
`rm -rf ./ds_rep1`;
`rm -rf ./trinity_out_dir`;
`rm -rf ./edgeR_isoforms`;
`rm -rf ./edgeR_genes`;
`rm -rf ./ds_rep1.stat`;
`rm -rf ./edgeR_trans`;
`rm -rf ./hs_rep1.stat`;
`rm -rf ./log_rep1.stat`;
`rm -rf ./plat_rep1.stat`;


exit(0);
