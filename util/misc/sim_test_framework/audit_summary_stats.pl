#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.audit_summary\n\n";

my $filename = $ARGV[0] or die $usage;

main: {

    my %yes_no_counter;
    my $total_reco = 0;
    my $total_ref = 0;
    my $total_FL = 0;

    open (my $fh, $filename) or die $!;
    while (<$fh>) {
        if (/ref_fa/) { next; }
        chomp;
        my @x = split(/\t/);
        my $yes_or_no = $x[4];
        
        $yes_no_counter{$yes_or_no}++;

        my $num_reco = $x[3];
        $total_reco += $num_reco;

        my $num_ref = $x[1];
        $total_ref += $num_ref;
        
        my $num_FL = $x[2];
        $total_FL += $num_FL;

    }
    close $fh;
    
    my $num_yes = $yes_no_counter{YES} || 0;
    my $num_no = $yes_no_counter{NO} ||0;

    my $num_additional = $total_reco - $total_FL;
    print join("\t", "#YES", "#NO", "#FL", "#TOT_Trans", "#Ref", "#additional") . "\n";
    print join("\t", $num_yes, $num_no, $total_FL, $total_reco, $total_ref, $num_additional) . "\n";
    
    exit(0);
}
