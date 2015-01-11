#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Cwd;

use Getopt::Long qw(:config no_ignore_case bundling);


$ENV{LC_ALL} = 'C';

my $usage = <<_EOUSAGE_;

################################################################################################
#
# Required:
#
# --genome      genome in fasta format
# --reads       reads in fasta format
#
# Optional:
#
# --blat_params  quote-delimited params to pass to blat, eg. "-q=rna -t=dna -maxIntron=1000"  (default: "-q=rna -t=dna")
#
# --top         number of top hits (default: 10)
# --min_per_ID  minimum percent identity (default: 95)
#
##############################################################################################


_EOUSAGE_

	;


my ($genome_fa, $reads_fa);

my $blat_params = "-q=rna -t=dna";
my $top_hits = 10;
my $min_per_ID = 95;



&GetOptions( 'genome=s' => \$genome_fa,
			 'reads=s' => \$reads_fa,
			 
			 'blat_params=s' => \$blat_params,
			 
			 'top=i' => \$top_hits,
			 'min_per_ID=i' => \$min_per_ID,
	);

unless ($genome_fa && $reads_fa) {
	die $usage;
}


{
	my @required_progs = qw (blat psl2sam.pl);

	foreach my $prog (@required_progs) {
		my $path = `which $prog`;
		unless ($path =~ /^\//) {
			die "Error, cannot locate required program: $prog";
		}
	}
}



main: {

	my $util_dir = "$FindBin::Bin/../util";
	
	my $cmd = "$util_dir/fasta_to_tab.pl $reads_fa > $reads_fa.tab";
	&process_cmd($cmd) unless (-s "$reads_fa.tab");

	$cmd = "sort -T . -S 2G -k 1,1 $reads_fa.tab > $reads_fa.tab.sort";
	&process_cmd($cmd) unless (-s "$reads_fa.tab.sort");

	# run blat
	$cmd = "blat $blat_params $genome_fa $reads_fa $reads_fa.psl";
	&process_cmd($cmd);

	# convert to sam
	$cmd = "psl2sam.pl -q 0 -r 0 $reads_fa.psl | sort -T . -S 2G -k 1,1 > $reads_fa.psl.sam";
	&process_cmd($cmd);


	# add the reads
	$cmd = "$util_dir/blat_sam_add_reads2.pl $reads_fa.psl.sam $reads_fa.tab.sort > $reads_fa.psl.sam.wReads";
	&process_cmd($cmd);

	## prune output to top matches:
	$cmd = "$util_dir/top_blat_sam_extractor.pl $reads_fa.psl.sam.wReads $top_hits $min_per_ID > $reads_fa.psl.sam.wReads.top";
	&process_cmd($cmd);

	$cmd = "$FindBin::Bin/cigar_tweaker $reads_fa.psl.sam.wReads.top $genome_fa > $reads_fa.psl.sam.wReads.top.tweaked";
	&process_cmd($cmd);

	$cmd = "sort -T . -S 2G -k 3,3 -k 4,4n $reads_fa.psl.sam.wReads.top.tweaked > $reads_fa.psl.sam.wReads.top.tweaked.coordSorted.sam";
	&process_cmd($cmd);

	
	
	exit(0);
}

####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}

	return;
}


