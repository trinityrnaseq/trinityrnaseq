#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
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

        
        if ($core_read_name ne $prev_read_name) {

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

    my %scaffold_to_reads;
    foreach my $read (@reads) {
        my $scaff_name = $read->get_scaffold_name();
        push (@{$scaffold_to_reads{$scaff_name}}, $read);
    }

    foreach my $scaff (keys  %scaffold_to_reads) {
        
        my @scaff_reads = @{$scaffold_to_reads{$scaff}};

        &process_scaffold_pairs(@scaff_reads);
    }


    return;
}




sub process_scaffold_pairs {
    my @reads = @_;
    
    my @left_reads;
    my @right_reads;

    foreach my $read (@reads) {
        
        my $read_name = $read->get_read_name();

        #print "processing: $read_name\n";

        if ($read_name =~ m|/1$|) {
            push (@left_reads, $read);
        }
        elsif ($read_name =~ m|/2$|) {
            push (@right_reads, $read);
        }
    }

    unless (@left_reads && @right_reads) {
        #print "\t** no pairs...\n";
        return;
    }
    
    my $left_read = shift @left_reads;
    my $left_aligned_pos = $left_read->get_aligned_position();
    my $scaffold_name = $left_read->get_scaffold_name();
        
    my $right_read = shift @right_reads;
    my $right_aligned_pos = $right_read->get_aligned_position();
    
    $left_read->set_mate_scaffold_name($scaffold_name);
    $left_read->set_mate_scaffold_position($right_aligned_pos);
    
    $right_read->set_mate_scaffold_name($scaffold_name);
    $right_read->set_mate_scaffold_position($left_aligned_pos);

    $left_read->set_paired(1);
    $right_read->set_paired(1);

    $left_read->set_proper_pair(1);
    $right_read->set_proper_pair(1);

    $left_read->set_first_in_pair(1);
    $right_read->set_second_in_pair(1);


    print $left_read->toString() . "\n";
    print $right_read->toString() . "\n";
    

    return;
}


    

        
    
