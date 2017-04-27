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
        
		my $transcribed_strand = $sam_entry->get_query_transcribed_strand($SS_lib_type);
        
		my $ofh = ($transcribed_strand eq '+') ? $plus_ofh : $minus_ofh;
		
		print $ofh $sam_entry->toString() . "\n";
		
	}
	
	close $plus_ofh;
	close $minus_ofh;

	exit(0);
}
