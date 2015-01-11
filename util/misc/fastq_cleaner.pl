#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 inputFile out.cleanReads out.malformedReads\n\n";

my $inputFile = $ARGV[0] or die $usage;
my $cleanReads = $ARGV[1] or die $usage;
my $malformedReads = $ARGV[2] or die $usage;


open (my $ofh_clean, ">$cleanReads") or die "Error, can't write to $cleanReads";
open (my $ofh_malformed, ">$malformedReads") or die "Error, cannot write to $malformedReads";

open (my $fh, $inputFile) or die "Error, cannot open $inputFile";

my $counter = 0;
my $num_clean = 0;
my $num_dirty = 0;

my @rec;

my $line = <$fh>;

while ($line) {

	if ($line =~ /^\@/) {
		$counter++;
		
		print STDERR "\r[$counter] [$num_clean clean] [$num_dirty dirty]       " if ($counter % 10000 == 0);
		
		push (@rec, $line);
		
		$line = <$fh>;
		for (1..3) {
			push (@rec, $line);
			$line = <$fh>;
		}
		
		my $record_text = join("", @rec);

		my $header = shift @rec;
		my $seq = shift @rec;
		my $qual_header = shift @rec;
		my $qual_line = shift @rec;
				
		chomp $header;
		chomp $seq if $seq;
		chomp $qual_header if $qual_header;
		chomp $qual_line if $qual_line;
		
		if ($header && $seq && $qual_header && $qual_line && 
			$qual_header =~ /^\+/ && length($seq) == length($qual_line)) {
			
			# can do some more checks here if needed to be sure that the lines are formatted as expected.
			print $ofh_clean join("\n", $header, $seq, $qual_header, $qual_line) . "\n";
			$num_clean++;
		}
		else {
			print $ofh_malformed $record_text;
			$num_dirty++;
		}
		@rec = ();
	} else {
		$line = <$fh>;
	}
	
}

exit(0);


		
		
		
		
