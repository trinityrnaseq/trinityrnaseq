#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 flag_number\n\n";

my $flag = $ARGV[0] or die $usage;


my @tokens;

my $pair_flag = ($flag & 0x0001) ? "PAIRED" : "unpaired";
push (@tokens, $pair_flag);
	
my $query_unmapped_flag = ($flag & 0x0004) ? "query_unmapped" : "QUERY_MAPPED";
push (@tokens, $query_unmapped_flag);

if ($query_unmapped_flag eq "QUERY_MAPPED") {
    my $query_strand = ($flag & 0x0010) ? "QUERY_REVERSE_STRAND" : "QUERY_FORWARD_STRAND";
    push (@tokens, $query_strand);
}
	
	
if ($pair_flag eq "PAIRED") {
	
    my $mapped_proper_pair_flag = ($flag & 0x0002) ? "MAPPED_PROPER_PAIR" : "not_mapped_proper_pair";
    push (@tokens, $mapped_proper_pair_flag);
    
    my $mate_mapped = ($flag & 0x0008) ? "mate_unmapped" : "MATE_MAPPED";
    push (@tokens, $mate_mapped);
    
    if ($mate_mapped eq "MATE_MAPPED") {
        my $mate_strand = ($flag & 0x0020) ? "MATE_REVERSE_STRAND" : "MATE_FORWARD_STRAND";
        push (@tokens, $mate_strand);
    }
		
    my $first_in_pair = ($flag & 0x0040) ? "FIRST_IN_PAIRx40" : "SECOND_IN_PAIRx40";
    push (@tokens, $first_in_pair);
		
    my $second_in_pair = ($flag & 0x0080) ? "SECOND_IN_PAIRx80" : "FIRST_IN_PAIRx80";
    push (@tokens, $second_in_pair);
}
	

my $primary_flag = ($flag & 0x0100) ? "notprimary" : "PRIMARY";
push (@tokens, $primary_flag);

my $fails_quality_checks = ($flag & 0x0200) ? "FAILED" : "";
push (@tokens, $fails_quality_checks) if $fails_quality_checks;

my $pcr_op_duplicate = ($flag & 0x0400) ? "PCROPDUP" : "";
push (@tokens, $pcr_op_duplicate) if $pcr_op_duplicate;

print "flag($flag) = " . join ("...", @tokens) . "\n";



exit(0);


