#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fastq_reader;
use File::Basename;

my $usage = "usage: $0 left.fq right.fq num_entries\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;
my $num_to_sample = $ARGV[2] or die $usage;

srand();

main: {

    ## do reservoir sampling on indices instead of the fastq records themselves, then just extract the selected ones.

    my $num_fq_entries = &get_num_fq_entries($left_fq);
    my @selected_entries = (1..$num_to_sample);
    
    for (my $i = $num_to_sample + 1; $i <= $num_fq_entries; $i++) {
        
        my $rand_pos = rand($i);
        if ($rand_pos < $num_to_sample) {
            $selected_entries[$rand_pos] = $i;
        }
    }

    my %indices_want = map { + $_ => 1 } @selected_entries;
    @selected_entries = (); # free

    
    my $left_fq_reader = new Fastq_reader($left_fq) or die("unable to open $left_fq for input");
    my $right_fq_reader = new Fastq_reader($right_fq) or die("unable to open $right_fq for input");;
    
    my $num_M_entries = $num_to_sample/1e6;
    $num_M_entries .= "M";
    my $base_left_fq = basename($left_fq);
    my $base_right_fq = basename($right_fq);
    
    open (my $left_ofh, ">$base_left_fq.$num_M_entries.fq") or die $!;
    open (my $right_ofh, ">$base_right_fq.$num_M_entries.fq") or die $!;
    
    my $fq_index = 0;
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
      
        $fq_index++;

        if ($indices_want{$fq_index}) {
            print $left_ofh $left_entry->get_fastq_record();
            print $right_ofh $right_entry->get_fastq_record();
        }
    }
    
    print STDERR " done.\n";
    
    close $left_ofh;
    close $right_ofh;
    
    exit(0);
}

####
sub get_num_fq_entries {
    my ($fastq_file) = @_;

    unless (-s $fastq_file) {
        die "Error, cannot locate file $fastq_file";
    }

    my $linecount;
    if ($fastq_file =~ /\.gz$/) {
        $linecount = `gunzip -c $fastq_file | wc -l`;
    }
    else {
        $linecount = `cat $fastq_file | wc -l`;
    }
    chomp $linecount;
    
    my $num_records = $linecount / 4;

    return($num_records);
}

