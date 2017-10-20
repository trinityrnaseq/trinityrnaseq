#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.fc\n";

my $fc_file = $ARGV[0] or die $usage;

open(my $fh, $fc_file) or die "Error, cannot open file $fc_file";
my $header_1 = <$fh>;
if ($header_1 !~ /^\#/) { die "Error, unexpected formatting of file: $fc_file";}
my $header_2 = <$fh>;
if ($header_2 !~ /^Geneid/) { die "Error, unexpected formatting of file: $fc_file";}
                          
while(<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $exon_num = $x[0];
    my $feature_count = $x[6];

    # round it off
    $feature_count = int($feature_count + 0.5);
    
    print "$exon_num\t$feature_count\n";
}
close $fh;

exit(0);


    
    
