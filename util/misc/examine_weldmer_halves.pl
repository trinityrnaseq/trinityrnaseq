#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 weldmer\n\n";

my $weldmer = $ARGV[0] or die $usage;

main: {
    
    if (-f $weldmer) {
        open (my $fh, $weldmer) or die $!;
        while (<$fh>) {
            my $weldmer = $_;
            chomp $weldmer;
            &check_weldmer($weldmer);
        }
        close $fh;
    }
    else {
        &check_weldmer($weldmer);
    }

    exit(0);
}


####
sub check_weldmer {
    my ($weldmer) = @_;
    
    my $weldmer_len = length($weldmer);
    
    my $left_substr = substr($weldmer, 0, int($weldmer_len/2));
    my $right_substr = substr($weldmer, int($weldmer_len/2));
    
    # &check_per_id($weldmer);
    &check_per_id($left_substr);
    &check_per_id($right_substr);
    
    
    return;
    
}


####
sub check_per_id {
    my ($weldmer) = @_;
    
    my $weldmer_len = length($weldmer);

    my $half_len = int($weldmer_len/2);
    #print "Half_len of " . length($weldmer_len) . " = $half_len\n";
    
    my @chars = split(//, $weldmer);

    my $max_ratio = 0;
   
    my $best_left = "";
    my $best_right = "";

    for (my $i = 0; $i < $half_len; $i++) {
        
        for (my $j = $i + 1; $j <= $half_len; $j++) {

            my $ref_pos = $i;
            my $other_pos = $j;
            
            my $count_same = 0;
            my $count_chars = 0;
            
            my $left_sub = "";
            my $right_sub = "";
            
            
            
            while ($other_pos <= $j + $half_len - 1) {
                
                $count_chars++;
                if ($chars[$ref_pos] eq $chars[$other_pos]) {
                    $count_same++;
                }
            
                $left_sub .= $chars[$ref_pos];
                $right_sub .= $chars[$other_pos];
                
                $ref_pos++;
                $other_pos++;
                
            }
            
            my $ratio = $count_same/$count_chars;
            if ($ratio > $max_ratio) {
                $max_ratio = $ratio;
                $best_left = $left_sub;
                $best_right = $right_sub;
            }
            
        }
    }
    
    $max_ratio = sprintf("%.3f", $max_ratio);
    print "$weldmer\t$best_left\t$best_right\t$max_ratio\n";
    
    return;
}
