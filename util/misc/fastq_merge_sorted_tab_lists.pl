#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: left.tab.sort right.tab.sort\n\n";

my $left_tab = $ARGV[0] or die $usage;
my $right_tab = $ARGV[1] or die $usage;


open (my $left_fh, $left_tab) or die $!;
open (my $right_fh, $right_tab) or die $!;

my $left_entry = <$left_fh>;
my $right_entry = <$right_fh>;

open (my $left_ofh, ">fixed.$$.left.fq") or die $!;
open (my $right_ofh, ">fixed.$$.right.fq") or die $!;

open (my $left_broken_ofh, ">broken.$$.left.fq") or die $!;
open (my $right_broken_ofh, ">broken.$$.right.fq") or die $!;

my $counter=0;

while (1) {
    
    $counter++;
    if ($counter % 1000000 == 0) {
        print STDERR "\r[$counter]   ";
    }
    
    unless ($left_entry && $right_entry) {
        last;
    }

    chomp $left_entry;
    chomp $right_entry;

    my ($left_acc, $left_seq, $left_qual) = split(/\t/, $left_entry);
    my ($right_acc, $right_seq, $right_qual) = split(/\t/, $right_entry);

    my $left_core = $left_acc;
    $left_core =~ s|/1$||;
    
    my $right_core = $right_acc;
    $right_core =~ s|/2$||;

    if ($left_core eq $right_core) {
        ## write entries
        
        if (length($left_qual) == length($left_seq)
            && 
            length($right_qual) == length($right_seq)) {
            
            ## AFAICT record looks good
            
            print $left_ofh join("\n", "\@$left_acc", $left_seq, "+", $left_qual) . "\n";
            print $right_ofh join("\n", "\@$right_acc", $right_seq, "+", $right_qual) . "\n";
        }
        else {
            print $left_broken_ofh $left_entry . "\n";
            print $right_broken_ofh $right_entry . "\n";
        }

        # reprime
        $left_entry = <$left_fh>;
        $right_entry = <$right_fh>;
    }
    else {
        if ($left_core lt $right_core) {
            print $left_broken_ofh $left_entry . "\n";
            $left_entry = <$left_fh>;
        }
        else {
            print $right_broken_ofh $right_entry . "\n";
            $right_entry = <$right_fh>;
        }
    }
}

close $left_ofh; 
close $right_ofh;
close $left_broken_ofh;
close $right_broken_ofh;


exit(0);


