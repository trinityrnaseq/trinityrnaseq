#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::RealBin or die "error, cannot cd to $FindBin::RealBin";




`rm -rf ./trinity_out_dir*`;
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
`rm -rf ./trinity_trim_and_norm_outdir`;

`rm -rf ./RSEM_outdir`;
`rm -rf ./eXpress_outdir`;
`rm -rf ./eXpress_PE`;
`rm -rf ./kallisto_PE`;
`rm -rf ./RSEM_PE`;
`rm -rf ./test_trinity_long_reads`;
`rm -rf ./test_trinity_bowtie2`;
`rm -rf ./trinity_test_samples_SE`;
`rm -rf ./trinity_test_samples_PE`;
`rm -rf ./trinity_piecemeal`;
`rm -rf ./trinity_test_samples_PE_min_kmer_cov_3`;
`rm -rf ./trinity_complete`;
`rm -rf ./trinity_incl_supertrans`;


exit(0);
