#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 fpkm.matrix pct_of_top_fpkm\n\n";

my $expr_matrix = $ARGV[0] or die $usage;
my $pct_of_top = $ARGV[1] or die $usage;

main: {

    open (my $fh, $expr_matrix) or die $!;
    
    my $header = <$fh>;
    chomp $header;
    $header =~ s/^\s+//;
    my @samples = split(/\t/, $header);
    while (<$fh>) {
        my ($gene, @vals) = split(/\t/);
        my $max_val = &get_max(@vals);
        
        my @within_range;
        for (my $i = 0; $i <= $#vals; $i++) {
            my $val = $vals[$i];
            my $pct_top = 100 - ($val/$max_val * 100);
            if ($pct_top <= $pct_of_top) {
                push (@within_range, $samples[$i]);
            }
        }
        my $sample_spec_type = join("|", sort @within_range);
        print "$gene\t$sample_spec_type\n";
    }
    

    close $fh;

    exit(0);
}

####
sub get_max {
    my @vals = @_;

    @vals = sort {$a<=>$b} @vals;

    my $max = pop @vals;

    return($max);
}
        
