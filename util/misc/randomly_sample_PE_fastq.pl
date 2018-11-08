#!/usr/bin/env perl

## util/fastQ_rand_subset.pl
## Purpose: Extracts out a specific number of random reads from an
##   input file, making sure that left and right ends of the selected
##   reads are paired
## Usage: $0 left.fq right.fq num_entries
##

## See:
## * https://en.wikipedia.org/wiki/Reservoir_sampling
##
## The current implementation requires the selected entries to be
## stored in memory, but passes over the input file(s) only once. If
## the array loading / reading is too slow or memory intensive, this
## could be implemented in a two-pass fashion by first getting
## indexes, and then getting the reads

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fastq_reader;
use File::Basename;
use List::Util qw(shuffle);
use Data::Dumper;

my $usage = "usage: $0 left.fq right.fq num_entries num_total_records\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;
my $num_entries = $ARGV[2] or die $usage;
my $num_total_records = $ARGV[3] or die $usage;


main: {

    srand();
    my @selected_indices = &get_random_indices($num_entries, $num_total_records);
    

    my %selected = map { + $_ => 1 } @selected_indices;

    @selected_indices = ();
    
    print STDERR "Selecting $num_entries entries...";
    
    my $left_fq_reader = new Fastq_reader($left_fq) or die("unable to open $left_fq for input");
    my $right_fq_reader = new Fastq_reader($right_fq) or die("unable to open $right_fq for input");;
    
    my $num_M_entries = $num_entries/1e6;
    $num_M_entries .= "M";
    my $base_left_fq = basename($left_fq);
    my $base_right_fq = basename($right_fq);

    
    open (my $left_ofh, ">$base_left_fq.$num_M_entries.fq") or die $!;
    open (my $right_ofh, ">$base_right_fq.$num_M_entries.fq") or die $!;
    
    my $counter = 0;
    while (my $left_entry = $left_fq_reader->next()) {
        my $right_entry = $right_fq_reader->next();
        
        unless ($left_entry && $right_entry) {
            die "Error, didn't retrieve both left and right entries from file ($left_entry, $right_entry) ";
        }
        unless ($left_entry->get_core_read_name() eq $right_entry->get_core_read_name()) {
            die "Error, core read names don't match: "
                . "Left: " . $left_entry->get_core_read_name() . "\n"
                . "Right: " . $right_entry->get_core_read_name() . "\n";
        }
        
        $counter++;

        if ($selected{$counter}) {
            
            print $left_ofh $left_entry->get_fastq_record();
            
            print $right_ofh $right_entry->get_fastq_record();
        
            delete($selected{$counter});

        }
    }

    print STDERR " done.\n";

    close $left_ofh;
    close $right_ofh;


    if (%selected) {
        die "Error, missing indices: " . Dumper(\%selected);
    }
        
    exit(0);
}



####
sub get_random_indices {
    my ($num_entries, $num_total_records) = @_;

    return( (shuffle 0..$num_total_records)[0..($num_entries-1)]);

}

