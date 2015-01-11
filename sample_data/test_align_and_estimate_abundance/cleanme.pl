#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 

                        reads.right.fq.gz
                        reads.left.fq.gz

                        run_eXpress.sh
                        run_rsem.sh

                        Trinity.fasta.gene_trans_map                        
                        Trinity.fasta.gz

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

`rm -rf ./eXpress_estimate`;
`rm -rf ./rsem_estimate`;


exit(0);
