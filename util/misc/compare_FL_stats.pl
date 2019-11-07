#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\tusage: $0 A.FL_selected B.FL_selected [RESTRICT_TO_FAILURES]\n\n";

my $A_stats_file = $ARGV[0] or die $usage;
my $B_stats_file = $ARGV[1] or die $usage;

my $RESTRICT_TO_FAILURES = $ARGV[2] || 0;

my $MIN_PER_LEN = 99;
my $MAX_PER_GAP = 1;

main: {
    
    my %A_stats = &parse_stats($A_stats_file);
    my %B_stats = &parse_stats($B_stats_file);

    my %trans = map { + $_ => 1 } (keys %A_stats, keys %B_stats);

    foreach my $trans_acc (keys %trans) {

        my $A_FL = $A_stats{$trans_acc} || ".";
        
        my $B_FL = $B_stats{$trans_acc} || ".";
        
        if ($RESTRICT_TO_FAILURES) {
            unless ($A_FL eq "." || $B_FL eq ".") {
                next;
            }
        }
        
        print join("\t", $A_FL, $B_FL) . "\n";
        
    }
    
    exit(0);
        
}

####
sub parse_stats {
    my ($stats_file) = @_;

    my %trans_to_stats;
    
    open (my $fh, $stats_file) or die $!;
    while (<$fh>) {
        chomp;

        my $line = $_;
        
        my @x = split(/\t/);
        

        my $trans = $x[0];
        my $reco = $x[1];

        $trans_to_stats{$trans} = $reco;

    }
    close $fh;

    return (%trans_to_stats);
}

