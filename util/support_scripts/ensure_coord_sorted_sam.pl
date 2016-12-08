#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\n\n\tusage: $0 coordinate_sorted.bam|sam\n\n";

my $bam_file = $ARGV[0] or die $usage;



my $num_records_to_validate = 1000;

main: {

    my $order_counter = 0;

    my $sam_reader =  new SAM_reader($bam_file);

       
    my $ordered_record_counter = 0;
    
    my $prev_record = $sam_reader->get_next();
    
    my %scaff_seen;
    $scaff_seen{ $prev_record->get_scaffold_name() } = 1;
    
    while (my $sam_entry = $sam_reader->get_next()) {
        
        
        if ($prev_record->get_scaffold_name() eq $sam_entry->get_scaffold_name()) {
            ## ensure coordinates are in order

            if ($sam_entry->get_scaffold_position() < $prev_record->get_scaffold_position()) {

                die "Error, read entries are out of order:\n" 
                    . $prev_record->get_original_line() . "\n"
                    . $sam_entry->get_original_line() . "\n"
                    . "\n\nBe sure to use a coordinate-sorted bam file\n";
            }
            elsif ($sam_entry->get_scaffold_position() > $prev_record->get_scaffold_position()) {
                # good, as we expect
                $ordered_record_counter++;
                
                if ($ordered_record_counter >= $num_records_to_validate) {
                    print STDERR "-appears to be a coordinate sorted bam file. ok.\n";
                    exit(0);
                }
            }
        }
        elsif ($scaff_seen{ $sam_entry->get_scaffold_name() } ) {
            die "Error, bam file doesn't appear to be coordinate sorted. Scaffold: " . $sam_entry->get_scaffold_name() . " is out of order. ";
        }
        
        $scaff_seen{ $sam_entry->get_scaffold_name() } = 1;
        
        $prev_record = $sam_entry;
    }
    
    
    die "ERROR: didn't find at least $num_records_to_validate BAM records properly ordered along a single scaffold. ";


}
    
    
        
    
