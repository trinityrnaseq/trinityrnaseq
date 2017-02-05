#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::RealBin or die "error, cannot cd to $FindBin::RealBin";



my @files_to_keep = qw (
cleanme.pl 
Makefile
align_and_estimate_tester_SINGLE_END.pl
samples.txt

test_Schizo.sh
schizo.samples.txt

test_Mouse.sh
mouse.samples.txt

test_Drosoph.sh
drosoph.samples.txt

);



my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (-f $file && ! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


`rm -rf ./RSEM-*`;
`rm -rf ./express-*`;
`rm -rf ./kallisto-*`;
`rm -rf ./salmon-*`;
`rm -rf Trinity.fasta.salmon*`;

`rm -rf mm9_ucsc_genomeStudio_genes.gtf*`;
`rm -rf SP2_GeneDB.290110.prot_coding.cds.nr*`;
`rm -rf ./Drosoph_coding.cdna.salmon_fmd.idx/`;
`rm -rf ./Drosoph_coding.cdna.salmon_quasi.idx`;

exit(0);
