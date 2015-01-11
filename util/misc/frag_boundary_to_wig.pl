#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.frag_coords (LEND|REND)\n\n";

my $file = $ARGV[0] or die $usage;
my $end = $ARGV[1] or die $usage;

unless ($end eq "LEND" || $end eq "REND") {
    die $usage;
}

main: {

    my @histogram;
    my $curr_scaffold = "";
    
    open (my $fh, $file) or die $usage;
    while (<$fh>) {
        chomp;
        my ($scaff, $read, $lend, $rend) = split(/\t/);
        if ($curr_scaffold ne $scaff) {
            &write_wig($curr_scaffold, \@histogram) if @histogram;
        
            # re-init
            $curr_scaffold = $scaff;
            @histogram = ();
        }
        
        my $pos = ($end eq "LEND") ? $lend : $rend;
        
        
        $histogram[$pos]++;
    }
    close $fh;
    
    if (@histogram) {
        &write_wig($curr_scaffold, \@histogram);
    }

    exit(0);
}

####
sub write_wig {
    my ($scaffold, $histogram_aref) = @_;

    print "variableStep chrom=$scaffold\n";

    for (my $i = 1; $i <= $#$histogram_aref; $i++) {
        my $val = $histogram_aref->[$i];
        if (defined $val) {
            print join("\t", $i, $val) . "\n";
        }
    }

    return;
}


