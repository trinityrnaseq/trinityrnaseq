#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 venn.txt\n\n";

my $venn = $ARGV[0] or die $usage;

my %class_to_att_count;

my %atts;

open (my $fh, $venn) or die "Error, cannot open file $venn";
while (<$fh>) {
	chomp;
	my ($feature, $class_list) = split(/\t/);

	my @classes = split(/,/, $class_list);
	
	my $att;
	
	my ($reference) = grep { /reference/ } @classes;

	my $count = scalar (@classes);
	if ($reference) {
		$count--;
		$att = "reference";
		@classes = grep { $_ !~ /reference/ } @classes;
	}
	
	
	elsif ( $count == 1) {
		$att = "unique";
	}
	elsif ($count > 1) {
		$att = "shared";
	}
	else {
		die "Error, have only reference...?";
	}
	
	$atts{$att}++;
	
	foreach my $class (@classes) {
		$class_to_att_count{$class}->{$att}++;
	}
}

close $fh;

my @sorted_atts = reverse sort keys %atts;

print "#class\t" . join("\t", @sorted_atts) . "\n";

my @classes = sort keys %class_to_att_count;


foreach my $class (@classes) {
	
	print $class;
	
	foreach my $att (@sorted_atts) {
		
		my $count = $class_to_att_count{$class}->{$att} || 0;
		print "\t$count";
	}
	
	print "\n";
}

exit(0);



		
	
	
