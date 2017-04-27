#!/usr/bin/env perl

## util/fastQ_rand_subset.pl
## Purpose: Extracts out a specific number of random reads from an
##   input file, making sure that left and right ends of the selected
##   reads are paired
## Usage: $0 left.fq right.fq num_entries
##
## This code implements reservoir sampling, which produces an unbiased
## selection regardless of the number of entries (assuming an
## appropriately uniform random distribution function).

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

my $usage = "usage: $0 single.fq num_entries\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $num_entries = $ARGV[1] or die $usage;

my @selected_left_entries = ();
my @selected_right_entries = ();

main: {

  print STDERR "Selecting $num_entries entries...";

  my $left_fq_reader = new Fastq_reader($left_fq) or die("unable to open $left_fq for input");
  
  my $num_M_entries = $num_entries/1e6;
  $num_M_entries .= "M";
  my $base_left_fq = basename($left_fq);

  open (my $left_ofh, ">$base_left_fq.$num_M_entries.fq") or die $!;

  srand();

  my $num_skipped = 0;
  my $num_output_entries = 0;
  my $num_entries_read = 0;

  my $left_entry = 0;
  my $right_entry = 0;

  while ($left_entry = $left_fq_reader->next()) {

      
    $num_entries_read++;

    if($num_entries_read <= $num_entries){
      # Populate reservoir with entries up to num_entries
      push(@selected_left_entries, $left_entry);
    
    } else {
      # Randomly replace elements in the reservoir
      # with a decreasing probability
      my $swapPos = int(rand($num_entries_read));

      if($swapPos < $num_entries){
        $selected_left_entries[$swapPos] = $left_entry;
      }
    }
  }

  if ($num_entries > $num_entries_read) {
    die "Error, num_entries $num_entries > total records available: $num_entries_read ";
  }

  # print selected entries to their respective files
  foreach my $entry (@selected_left_entries){
    print $left_ofh $entry->get_fastq_record();
  }
  
  print STDERR " done.\n";

  close $left_ofh;
  
  exit(0);
}


