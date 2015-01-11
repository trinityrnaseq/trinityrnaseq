#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.sam [include_SAM=0]\n\n";

my $SAM_FILE = $ARGV[0] or die $usage;
my $INCLUDE_SAM = $ARGV[1];


my $fh;
if ($SAM_FILE eq "-") {
	$fh = \*STDIN;
}
else {
	if ($SAM_FILE =~ /\.bam$/) {
        open ($fh, "samtools view $SAM_FILE |") or die "Error, $!";
    }
    else {
        open ($fh, $SAM_FILE) or die "Error, cannot open file $SAM_FILE";
    }
}


while (<$fh>) {
	chomp;

	my $line = $_;
	
	if (/^\@/) { next; } # skip header
	my @x = split(/\t/);
	
	my $read_acc = $x[0];
	my $flag = $x[1];
	my $qual_string = $x[10];
	my $target = $x[2];
    
	my @tokens;
	
	my $pair_flag = ($flag & 0x0001) ? "PAIRED" : "notpaired";
	push (@tokens, $pair_flag);
	
	my $query_unmapped_flag = ($flag & 0x0004) ? "qun" : "QM";
	push (@tokens, $query_unmapped_flag);
	
	if ($query_unmapped_flag eq "QM") {
		my $query_strand = ($flag & 0x0010) ? "QREV" : "QFWD";
		push (@tokens, $query_strand);
	}
	
	
	if ($pair_flag eq "PAIRED") {
	
		my $mapped_proper_pair_flag = ($flag & 0x0002) ? "MAPPEDPROPERPAIR" : "notmappedproperpair";
		push (@tokens, $mapped_proper_pair_flag);
		
		my $mate_mapped = ($flag & 0x0008) ? "mateunmapped" : "MATEMAPPED";
		push (@tokens, $mate_mapped);
		
		if ($mate_mapped eq "MM") {
			my $mate_strand = ($flag & 0x0020) ? "MREV" : "MFWD";
			push (@tokens, $mate_strand);
		}
		
		my $first_in_pair = ($flag & 0x0040) ? "FIRSTx40" : "SECONDx40";
		push (@tokens, $first_in_pair);
		
		my $second_in_pair = ($flag & 0x0080) ? "SECONDx80" : "FIRSTx80";
		push (@tokens, $second_in_pair);
	}
	
	
	my $primary_flag = ($flag & 0x0100) ? "notprimary" : "PRIMARY";
	push (@tokens, $primary_flag);
	
	my $fails_quality_checks = ($flag & 0x0200) ? "FAILED" : "";
	push (@tokens, $fails_quality_checks) if $fails_quality_checks;
	
	my $pcr_op_duplicate = ($flag & 0x0400) ? "PCROPDUP" : "";
	push (@tokens, $pcr_op_duplicate) if $pcr_op_duplicate;
	
	if ($INCLUDE_SAM) {
		print $line;
	}
	else {
		print "$read_acc\t$target";
	}
	print "\t" . join ("_", @tokens) . "\n";
	
}

exit(0);


