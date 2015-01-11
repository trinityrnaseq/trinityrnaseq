#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (
                        Trinity.seq_lengths
                        Trinotate_report.xls
                        Trinotate_report.xls.trans.gene_ontology
                        cleanme.pl
                        ds_induced_vs_log.factors
                        hs_induced_vs_log.factors
                        runMe.sh
                        );


my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}



exit(0);
