#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: all_diff_expr_results.txt min_abs_log_FC max_FDR\n\n";

my $dat_file = $ARGV[0] or die $usage;
my $min_logFC = $ARGV[1] or die $usage;
my $max_FDR = $ARGV[2] or die $usage;


my %best;

open (my $fh, $dat_file) or die $!;
while (<$fh>) {
    my $line = $_;
    chomp;
    if (/^\#/) { next; }
    my @x = split(/\t/);
    
    my $acc = $x[2];
    
    my $logfc = $x[4];
    my $FDR = $x[6];
    
    if (abs($logfc) >= $min_logFC && $FDR <= $max_FDR) {
        
        if (my $best_entry = $best{$acc}) {
            if ($best_entry->{FDR} > $FDR) {
                
                # replace it
                $best{$acc} = { line => $line,
                                FDR => $FDR,
                };
            }
        }
        else {
            # init
                $best{$acc} = { line => $line,
                                FDR => $FDR,
                };
        }
    }
}


foreach my $entry (values %best) {
    print $entry->{line};
}


exit(0);


