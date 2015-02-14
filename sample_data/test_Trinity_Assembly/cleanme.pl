#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (cleanme.pl 
                        README
                        reads.left.fq.gz
                        reads.right.fq.gz
                        runMe.sh
                        __indiv_ex_sample_derived
                        misc_run_tests
                        run_abundance_estimation_procedure.sh

                        reads2.left.fq.gz
                        reads2.right.fq.gz

                        longReads.fa
                        Makefile
 test_FL.sh

                        __test_runMe_with_jaccard_clip.sh
                        __test_runMe_with_SE_reads.sh

   runRecTrinity.sh
                        );


my %keep = map { + $_ => 1 } @files_to_keep;


`rm -rf ./trinity_out_dir` if (-d "trinity_out_dir");
`rm -rf ./bowtie_out` if (-d "bowtie_out");
`rm -rf ./RSEM.stat` if (-d "RSEM.stat");
`rm -rf ./trinity_single_outdir` if (-d "trinity_single_outdir");
`rm -rf ./bowtie_single_outdir` if (-d "bowtie_single_outdir");

`rm -rf ./normalized_reads_test` if (-d "normalized_reads_test");
`rm -rf ./collectl` if (-d "collectl");
`rm -rf ./trin_w_jaccard` if (-d "trin_w_jaccard");
`rm -rf ./trin_SE_outdir` if (-d "trin_SE_outdir");
`rm -rf ./trin_test_no_qtrim`;

`rm -rf ./trinity_test_no_qtrim` if (-d "trinity_test_no_qtrim");
`rm -rf ./trinity_w_jaccard` if (-d "trinity_w_jaccard");

`rm -rf ./__test_trinity*`;


foreach my $file (<*>) {
	
	if (! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
