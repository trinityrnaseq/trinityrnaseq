#!/usr/bin/env perl

use strict;
use warnings;


## Given a subset of the venn, extracts all entries from the original files containing the pairwise DE data.

my $usage = "\n\nusage: $0 venn.selected pw1.summary pw2.summary ...\n\n";

my $selected_file = $ARGV[0] or die $usage;
shift @ARGV;

my @pw_summaries = @ARGV;
unless (scalar @pw_summaries > 1) { die $usage; }


main: {

    my %keys_want;
    {
        open (my $fh, $selected_file) or die "Error, cannot open file $selected_file";
        while (<$fh>) {
            chomp;
            my ($acc, $sampleA, $sampleB, @rest) = split(/\t/);
            
            my $key = join("$;", $acc, sort ($sampleA, $sampleB));

            $keys_want{$key} = 1;
            
        }
        close $fh;
    }

    foreach my $pw_summary (@pw_summaries) {
        
        open (my $fh, $pw_summary) or die $!;
        while (<$fh>) {
            my  $line = $_;
            chomp;
            
            my ($acc, $sampleA, $sampleB, @rest) = split(/\t/);
            my $key = join("$;", $acc, sort ($sampleA, $sampleB));

            
            if ($keys_want{$key}) {
                print $line;
            }
        }
        close $fh;

    }

    exit(0);
}


