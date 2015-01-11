#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 file.sam SS_lib_type={F,R,FR,RF}\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $SS_lib_type = $ARGV[1] or die $usage;

unless ($SS_lib_type =~ /^(F|R|FR|RF)$/) {
	die "Error, SS_lib_type must be F, R, FR, or RF";
}


main: {

	my $sam_reader = new SAM_reader($sam_file);


	my $num_verified = 0;
	my $num_conflicted = 0;

	my $counter = 0;

	while ($sam_reader->has_next()) {
		
		$counter++;
		print STDERR "\r[$counter]  " if $counter % 10000 == 0;

		#if ($counter % 300000 == 0) { last; } ## debugging
		
		my $sam_entry = $sam_reader->get_next();

		my $aligned_strand = $sam_entry->get_query_strand();
		my $opposite_strand = ($aligned_strand eq '+') ? '-' : '+';
		
		my $transcribed_strand = $aligned_strand;

		if ($sam_entry->is_paired()) {
			
			unless ($SS_lib_type =~ /^(FR|RF)$/) {
				die "Error, SS_lib_type: $SS_lib_type is not compatible with paired reads";
			}
			

			if ($sam_entry->is_first_in_pair()) {
				
				if ($SS_lib_type eq "RF") {
					$transcribed_strand = $opposite_strand;
				}
			}

			else {
				# second in pair
				if ($SS_lib_type eq "FR") {
					$transcribed_strand = $opposite_strand;
				}
			}
		}
		else {
			## Unpaired or Single Reads
			if ($SS_lib_type eq "R") {
				$transcribed_strand = $opposite_strand;
			}
		}

		my $sam_text = $sam_entry->toString();

		my @x = split(/\t/, $sam_text);

		my $strand_flag = "XS:A:$transcribed_strand";

		if (my ($XS_A) = grep { $_ =~ /^XS:A:/ } @x) {
			if ($XS_A eq $strand_flag) {
				$num_verified++;
			}
			else {
				$num_conflicted++;
				print STDERR "Conflicting strand info: existing = $XS_A, expected = $strand_flag\n";
				@x = grep { $_ !~ /^XS:A:/ } @x;
				push (@x, $strand_flag);
			}
			
		}
		else {
			push (@x, $strand_flag);
		}
		
		
		print join("\t", @x) . "\n";
	}

	if ($num_conflicted || $num_verified) {
		print STDERR "Percent entries with conflicted transcribed strand assignments = " 
			. sprintf("%.2f", $num_conflicted / ($num_conflicted + $num_verified) * 100) . "\n";
	}

	exit(0);

}
