#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Basename;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

##################################################################
#
# Required:
#
#     --target <string>    : transcript sequences (eg. 'Trinity.fasta')
#
#  If paired reads:
#
#     --left  <string>    :left reads
#     --right <string>    :right reads
# 
#  Or, if unpaired reads:
#
#      --single <string>   :single reads
#
#  
#  --paired_fragment_length <int>  :size of a read pair insert (def=300)
#
#
#############################################################################################################



_EOUSAGE_

	;



my ($target_fasta, $left_file, $right_file, $single_file, $SS_lib_type, $paired_fragment_length);

# defaults:
$paired_fragment_length = 300;


&GetOptions( 
			 
			 ## general opts
			 "target_fasta=s" => \$target_fasta,
			 "left=s" => \$left_file,
			 "right=s" => \$right_file,
			 "single=s" => \$single_file,
			 
			 "SS_lib_type=s" => \$SS_lib_type,
			 
			 "paired_fragment_length=i" => \$paired_fragment_length,
			 );


## Check options set:

unless  ( ($left_file && $right_file) || $single_file) {
	die $usage;
}


main: {

	my $cmd = "ln -s $target_fasta TARGET.fa";
	&process_cmd($cmd) unless (-e "TARGET.fa");
	
	$cmd = "bowtie-build TARGET.fa TARGET";
	&process_cmd($cmd) unless (-s "TARGET.1.ebwt");
	
	if ($left_file && $right_file) {
		$cmd = "tophat --bowtie1 -i 5 -r $paired_fragment_length TARGET $left_file $right_file";
		&process_cmd($cmd);
	}
	else {
		$cmd = "tophat --bowtie1 -i 5 TARGET $single_file";
		&process_cmd($cmd);
	}


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
