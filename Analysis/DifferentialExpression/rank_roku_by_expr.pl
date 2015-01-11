#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 roku.dat expr.dat\n\n";

my $roku_dat = $ARGV[0] or die $usage;
my $expr_dat = $ARGV[1] or die $usage;


my %acc_to_info;

main: {

    {
        open (my $fh, $roku_dat) or die $!;
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my ($acc, $score) = @x;
            $acc_to_info{$acc}->{roku} = $score;
            $acc_to_info{$acc}->{acc} = $acc;
        }
    }

    {

        open (my $fh, $expr_dat) or die $!;
        my $header = <$fh>;
        while (<$fh>) {
            chomp;
            my ($acc, $expr) = split(/\t/);
            $acc_to_info{$acc}->{expr} = $expr;
        }
    }

    my %bins;

    foreach my $entry (values %acc_to_info) {

        my $bin = int($entry->{roku}/0.5);

        push (@{$bins{$bin}}, $entry);
    }

    foreach my $bin (sort {$a<=>$b} keys %bins) {

        my @entries = @{$bins{$bin}};

        @entries = reverse sort {$a->{expr}<=>$b->{expr}} @entries;

        
        foreach my $entry (@entries) {
            print "$bin\t$entry->{acc}\t$entry->{expr}\t$entry->{roku}\n";
        }
        print "\n";
    }


    exit(0);
}


