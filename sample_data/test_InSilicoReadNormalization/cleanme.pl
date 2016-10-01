#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::RealBin or die "error, cannot cd to $FindBin::RealBin";



my @files_to_keep = qw (cleanme.pl 
                        test_PE_disordered.sh
                        test_PE_perfect_ordering.sh
                        test_PE_perfect_ordering.w_base_cov_stats.sh
                        test_SE_normalization.sh
                        Makefile
                        
                        test_PE_normalization.w_base_cov_stats.sh
                        test_PE_normalization.sh
test_PE_normalization.mult_read_sets.sh
test_PE_normalization_useFASTA.sh
                        
                        );

my %keep = map { + $_ => 1 } @files_to_keep;


`rm -rf ./tmp_normalized_reads/`;
`rm -rf ./single_tmp_norm_reads/`;
`rm -rf ./tmp_PE_norm_dir`;
`rm -rf ./test_multi_read_sets_norm_outdir/`;
`rm -rf ./tmp_PE_norm_FA_dir`;

foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
