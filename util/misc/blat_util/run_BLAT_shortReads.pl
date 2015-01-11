#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 genome_db transcript_db maxIntron outFile\n\n";

my $genome_db = $ARGV[0] or die $usage;
my $transcript_db = $ARGV[1] or die $usage;
my $max_intron = $ARGV[2] or die $usage;
my $outFile = $ARGV[3] or die $usage;

main: {

	my $ooc_cmd = "blat -t=dna -q=rna -maxIntron=$max_intron -makeOoc=11.ooc $genome_db $transcript_db $outFile";
	unless (-s "11.ooc") {
		&process_cmd($ooc_cmd);
	}


	my $blat_cmd = "blat -t=dna -q=rna -maxIntron=$max_intron -ooc=11.ooc $genome_db $transcript_db $outFile";
	&process_cmd($blat_cmd);
	
	exit(0);
}

sub process_cmd {
	my ($cmd) = @_;
	
	my $ret = system($cmd);
	if ($ret) {
		die "Error, $cmd died with ret $ret";
	}
	
	return;
}
