#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 min_len min_per_len min_per_id file.per_len_and_id ...\n";



my $min_len = shift @ARGV;
my $min_per_len = shift @ARGV;
my $min_per_id = shift @ARGV;
my @files = @ARGV;

unless ($min_len && $min_per_len && $min_per_id && @files) {
    die $usage;
}

foreach my $file (@files) {

    my $total_considered = 0;
    my $total_OK = 0;

    open (my $fh, $file) or die "Error, cannot open file $file";
    while (<$fh>) {
        chomp;
        my ($acc, $match_len, $trans_len, $per_len, $per_id) = split(/\t/);

        if ($match_len >= $min_len) {
            
            $total_considered++;
            
            if ($per_len >= $min_per_len && $per_id >= $min_per_id) {
                
                $total_OK++;
            }
        }

    }
    my $percent_OK = sprintf("%.2f", $total_OK / $total_considered);
    
    print "$file\t$total_considered\t$total_OK\t$percent_OK\n";
}

exit(0);

