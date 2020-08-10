#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use SAM_reader;
use SAM_entry;
use File::Basename;

my $usage = "usage: $0 alignments.sam SS_lib_type=F,R,FR,RF\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $SS_lib_type = $ARGV[1] or die $usage;

unless ($SS_lib_type =~ /^(F|R|FR|RF)$/) {
	die $usage;
}
	


main: {
	
    my $sam_reader = new SAM_reader($sam_file);

    my $sam_basename = basename($sam_file);
    
    
    open (my $plus_ofh, ">$sam_basename.+.sam") or die "Error, cannot write to $sam_basename.+.sam";
	open (my $minus_ofh, ">$sam_basename.-.sam") or die "Error, cannot write to $sam_basename.-.sam";
	
	while (my $sam_entry = $sam_reader->get_next()) {

        if ($sam_entry->is_query_unmapped()) {
            next;
        }
        
		my $transcribed_strand = &get_transcribed_strand($sam_entry);
        
		my $ofh = ($transcribed_strand eq '+') ? $plus_ofh : $minus_ofh;
		
		print $ofh $sam_entry->toString() . "\n";
		
	}
	
	close $plus_ofh;
	close $minus_ofh;

	exit(0);
}




####
sub get_transcribed_strand {
    my ($sam_entry) = @_;
    
    unless ($SS_lib_type) {
        confess "Error, SS_lib_type required as a parameter, possible values: RF,FR,F,R " . $sam_entry->toString();
    }
        
    my $aligned_strand = $sam_entry->get_query_strand();
    my $opposite_strand = ($aligned_strand eq '+') ? '-' : '+';
    
    my $transcribed_strand;
    
    if (! $sam_entry->is_paired()) {
        
        my $is_long_read = &get_long_read_status($sam_entry);
        
        if ($is_long_read) {
            return($aligned_strand); # should already be oriented in the forward orientation with respect to input sequence.
        }
        
        ## UNPAIRED or SINGLE READS
        unless ($SS_lib_type =~ /^(F|R)$/) {
            confess "Error, cannot have $SS_lib_type library type with unpaired reads " . $sam_entry->toString();
        }
        
        $transcribed_strand = ($SS_lib_type eq "F") ? $aligned_strand : $opposite_strand;
    }
    else {
        
        
        ## paired RNA-Seq reads:  left fragment is on the 3' end revcomped, and right fragment is at the 5' end sense strand.
        
                    
        unless ($SS_lib_type =~ /^(FR|RF)$/) {
            confess "Error, cannot have $SS_lib_type library type with paired reads " . $sam_entry->toString();
        }
                
        if ($sam_entry->is_first_in_pair()) {
            $transcribed_strand = ($SS_lib_type eq "FR") ? $aligned_strand : $opposite_strand;
        }
        else {
            # second pair
            $transcribed_strand = ($SS_lib_type eq "FR") ? $opposite_strand : $aligned_strand;
        }
    }
    
    
    return($transcribed_strand);
    
}



####
sub get_long_read_status {
    my ($sam_entry) = @_;
    

    my @fields = $sam_entry->get_fields();
    @fields = @fields[11..$#fields];
    
    if (grep { /^RG:Z:PBLR$/ } @fields) {
        return(1); 
    }
    else {
        return(0);
    }
}
