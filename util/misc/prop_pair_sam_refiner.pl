#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\n\nusage: name_sorted_paired_reads.sam\n\n";

my $sam_file = $ARGV[0] or die $usage;

main: {


    my $prev_read_name = "";
    my $prev_scaff_name = "";

    my @reads;

    my $sam_reader = new SAM_reader($sam_file);
    while ($sam_reader->has_next()) {
        

        my $read = $sam_reader->get_next();
        
        my $scaff_name = $read->get_scaffold_name();
        my $core_read_name = $read->get_core_read_name();
        
        if ($scaff_name ne $prev_scaff_name ||  $core_read_name ne $prev_read_name) {

            if (@reads) {
                &process_pairs(@reads);
                @reads = ();
            }
        }
        
        push (@reads, $read);
        

        $prev_read_name = $core_read_name;
        $prev_scaff_name = $scaff_name;
        
        
    }
    
    &process_pairs(@reads);


    exit(0);
}

####
sub process_pairs {
    my (@reads) = @_;


    my @left_reads;
    my @right_reads;

    foreach my $read (@reads) {
        if ($read->is_first_in_pair()) {
            push (@left_reads, $read);
        }
        elsif ($read->is_second_in_pair()) {
            push (@right_reads, $read);
        }
    }

    unless (@left_reads && @right_reads) {
        die "Error, dont have pairs!";
    }
    
    my $left_read = shift @left_reads;
    my $aligned_pos = $left_read->get_aligned_position();

    my $right_read = undef;
    while ((! defined($right_read)) && @right_reads) {
        my $read = shift @right_reads;
        if ($read->get_mate_scaffold_position() == $aligned_pos) {
            $right_read = $read;
            last;
        }
    }

    unless ($left_read && $right_read) {
        die "Error, couldn't match pairs.";
    }

    print $left_read->toString() . "\n";
    print $right_read->toString() . "\n";
    

    return;
}


    

        
    
