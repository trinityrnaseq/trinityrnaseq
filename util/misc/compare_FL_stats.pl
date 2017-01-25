#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\tusage: $0 before.stats after.stats\n\n";

my $before_stats_file = $ARGV[0] or die $usage;
my $after_stats_file = $ARGV[1] or die $usage;

my $MIN_PER_LEN = 99;
my $MAX_PER_GAP = 1;

main: {
    
    my %before_stats = &parse_stats($before_stats_file);
    my %after_stats = &parse_stats($after_stats_file);

    my %trans = map { + $_ => 1 } (keys %before_stats, keys %after_stats);

    foreach my $trans_acc (keys %trans) {

        print "// $trans_acc\n";

        my $before = $before_stats{$trans_acc} || "";
        print "# before\n$before\n";

        my $after = $after_stats{$trans_acc} || "";
        print "# after\n$after\n\n";

    }

    exit(0);
        
}

####
sub parse_stats {
    my ($stats_file) = @_;

    my %tran_to_stats;
    
    open (my $fh, $stats_file) or die $!;
    while (<$fh>) {
        chomp;

        my $line = $_;
        
        my @x = split(/\t/);
        

        my $acc = $x[0];
        my $per_gap = $x[3];
        my $per_len = $x[10];

        if ($per_len >= $MIN_PER_LEN && $per_gap <= $MAX_PER_GAP) {

            $line .= "\t**";
        }

        $tran_to_stats{$acc} .= "$line\n";

    }
    close $fh;
    
    return(%tran_to_stats);
    
}


