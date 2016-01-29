#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $usage = "usage: $0 reads.list [DS_flag=0]\n\n";

my $reads_file = $ARGV[0] or die $usage;
my $DS_flag = $ARGV[1];

open (my $fh, $reads_file) or die "Error, cannot open file $reads_file";
while (<$fh>) {
	my $file = $_;
	chomp $file;

	my $iworm_outfile = "$file.iworm";
	
	my $cmd = "$FindBin::RealBin/../Inchworm/bin/inchworm --reads $file --run_inchworm ";
	
	if ($DS_flag) {
		$cmd .= "--DS ";
	}
	
	$cmd .= "> $iworm_outfile";
	
	print "$cmd\n";
}

exit(0);




		
