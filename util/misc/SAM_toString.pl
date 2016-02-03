#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: $0 file.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;


main: {
    my $sam_reader = new SAM_reader($sam_file);
    
    while (my $sam_entry = $sam_reader->get_next()) {
        my ($genome_coords_aref, $read_coords_aref) = $sam_entry->get_alignment_coords();
        
        my $read_name = $sam_entry->get_read_name();
        my $scaffold = $sam_entry->get_scaffold_name();
        
        print join("\t", $scaffold, &get_coord_string($genome_coords_aref),
                   $read_name, &get_coord_string($read_coords_aref)) . "\n";
    }
    
    exit(0);
}


####
sub get_coord_string {
    my ($coords_aref) = @_;

    my $text = "";
    
    foreach my $coordset (@$coords_aref) {
        my ($lend, $rend) = @$coordset;
        
        if ($text) {
            $text .= ",";
        }
        
        $text .= "$lend-$rend";
    }

    return($text);
}
