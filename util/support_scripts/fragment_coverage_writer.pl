#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 file.sam.frag_coords\n\n";

my $frag_coords_file = $ARGV[0] or die $usage;

main: {
	
	my %scaffold_to_coverage; # will retain all coverage information.
	
    
    open (my $fh, $frag_coords_file) or die "Error, cannot open file $frag_coords_file";
	
	my $current_scaff = undef;
	
	my $counter = 0;
	while (my $line = <$fh>) {

        chomp $line;
        
		$counter++;

        my ($scaff, $frag_name, $lend, $rend) = split(/\t/, $line);

		if ($counter % 1000 == 0) {
			print STDERR "\r[$counter lines read]  scaff:$scaff lend:$lend  rend:$rend     ";
		}
        		
		if (%scaffold_to_coverage && $scaff ne $current_scaff) {
			&report_coverage(\%scaffold_to_coverage);
			%scaffold_to_coverage = ();
		}
		
		$current_scaff = $scaff;
        
		&add_coverage($current_scaff, [$lend, $rend], \%scaffold_to_coverage);
		

	}

	if (%scaffold_to_coverage) {
		&report_coverage(\%scaffold_to_coverage);
	}
		

	exit(0);

}


####
sub report_coverage {
	my ($scaffold_to_coverage_href) = @_;
	
	## output the coverage information:
	foreach my $scaffold (sort keys %$scaffold_to_coverage_href) {
		
        unless (defined($scaffold) && $scaffold =~ /\w/) { next; }
        
		print "variableStep chrom=$scaffold\n";
		
		my $coverage_aref = $scaffold_to_coverage_href->{$scaffold};

        my $first_cov_pos = &find_first_cov_pos($coverage_aref);
		
		for (my $i = $first_cov_pos; $i <= $#$coverage_aref; $i++) {
			my $cov = $coverage_aref->[$i] || 0;
			
			print "$i\t$cov\n";
		}
		
	}
	
	return;
}


####
sub add_coverage {
	my ($scaffold, $coords_aref, $scaffold_to_coverage_href) = @_;
	
    my ($lend, $rend) = @$coords_aref;
	
    ## add coverage:
    for (my $i = $lend; $i <= $rend; $i++) {
        $scaffold_to_coverage_href->{$scaffold}->[$i]++;
    }
    

	return;
}
	
####
sub find_first_cov_pos {
    my ($coverage_aref) = @_;

    for (my $i = 0; $i <= $#$coverage_aref; $i++) {
        if ($coverage_aref->[$i]) {
            return($i);
        }
    }
    
    return($#$coverage_aref);
}


