#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 
                        runMe.sh
                        samples_n_reads_decribed.txt
                        
                        Makefile
                        );


my %keep = map { + $_ => 1 } @files_to_keep;

`rm -rf edgeR_\*`;
#`find ./trinity_out_dir -type f -delete`;
#`rm -rf ./trinity_out_dir`;
`rm -f rnaseq_reads/*fq`;
`rm -rf collectl/`;
`rm -rf *.stat/`;
`rm -rf ./read_content_analysis/`;

foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
