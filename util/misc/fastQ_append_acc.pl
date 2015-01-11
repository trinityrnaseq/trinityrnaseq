#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.fastq append_no\n\n";

my $fastq_file = $ARGV[0] or die $usage;
my $append_no = $ARGV[1] or die $usage;


my @lines;

my $counter = 0;
open (my $fh, $fastq_file) or die "Error, cannot open file $fastq_file\n";
while (<$fh>) {
    chomp;
    $counter++;
    push (@lines, $_);
    
    if ($counter % 4 == 0) {
        if ($lines[0] !~ /^\@/) {
            die "Error, fastq record doesn't start with a header line as expected: " . join("\n", @lines);
        }

        my @x = split(/\s+/, $lines[0]);
        $x[0] .= "/$append_no";
        $lines[0] = join(" ", @x);
        print join("\n", @lines) . "\n";
        @lines = ();
    }

}
close $fh;

exit(0);

