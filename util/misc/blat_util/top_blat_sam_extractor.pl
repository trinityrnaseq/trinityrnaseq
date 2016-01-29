#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../../PerlLib");

use SAM_reader;
use SAM_entry;


my $usage = "usage: $0 blat.nameSorted.sam [num_top_hits=20] [min_per_ID=0]\n\n";

# note min perID is based on read length and blat2sam.pl default scoring with -q 0 -r 0

my $blat_sam = $ARGV[0] or die $usage;
my $num_top_hits = $ARGV[1] || 20;
my $min_per_ID = $ARGV[2] || 0;

main: {
	
	my $sam_reader = new SAM_reader($blat_sam);
	
	while ($sam_reader->has_next()) {
		
		my @entries;
		
		my $sam_entry = $sam_reader->get_next();
		push (@entries, $sam_entry);
		
		while ($sam_reader->has_next() 
			   &&
			   $sam_reader->preview_next()->get_read_name() eq $sam_entry->get_read_name()) {

			push (@entries, $sam_reader->get_next());
		}

		&report_top_hits(@entries);
	}

	exit(0);
	
}


####
sub report_top_hits {
	my @entries = @_;

    @entries = &get_top_hits(@entries);
	

	foreach my $entry (@entries) {
		print $entry->toString() . "\n";
	}
	
	return;
}



####
sub get_top_hits {
	my @entries = @_;

	my @structs;

	foreach my $entry (@entries) {
		my @fields = $entry->get_fields();

		my $seq_length = length($entry->get_sequence());
		
		my $min_score = $seq_length -  ( ( 1 - ($min_per_ID / 100)) * $seq_length * 3);  # perfect_score - (num_mismatches * mismatch_penalty)
				
		my $alignment_field = $fields[11];
		
		$alignment_field =~ /AS:i:(\d+)/ or die "Error, no score reported for " . $entry->toString();
		
		my $score = $1;
		unless (defined $score) {
			die "Error, no score for " . $entry->toString();
		}

        #print "MIN_score: $min_score vs. score: $score\n";
		
		if ($score >= $min_score) {
			
			push (@structs, { entry => $entry,
							  score => $score,
					  }
				);
			
		}
	}
	
	@structs = reverse sort {$a->{score}<=>$b->{score}} @structs;

    
	
    if ($num_top_hits < 0 && scalar(@structs) > 1) {
        # multiply mapped reads
        # ignoring entry.
        return();
    }
    
	elsif ($num_top_hits > 0 && scalar (@structs) > $num_top_hits) {
		@structs = @structs[0..$num_top_hits-1];
	}
	
	## unwrap:
	my @ret;
	if (@structs) {
	   my $top_struct = shift @structs;
		
		my $top_score = $top_struct->{score};
	   
	   push (@ret, $top_struct->{entry});
	   
		foreach my $struct (@structs) {
			if ($struct->{score} == $top_score) {
				push (@ret, $struct->{entry});
			}
			else {
				last;
			}
		}
    }
	
	return(@ret);
}
			 
