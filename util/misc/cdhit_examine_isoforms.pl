#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\tusage: $0 cd-hit.clstr\n\n";

my $cdhit_file = $ARGV[0] or die $usage;

main: {

    my $num_bad_clusters = 0;
    
    my $cluster;
    my @trans;
    
    open(my $fh, $cdhit_file) or die $!;
    while (<$fh>) {
        chomp;
        if (/^>/) {
            if (@trans) {
                $num_bad_clusters += &examine_cluster($cluster, \@trans);
            }
            $cluster = $_;
            @trans = ();
        }
        else {
            push (@trans, $_);
        }
    }
    close $fh;

    if (@trans) {
        $num_bad_clusters += &examine_cluster($cluster, \@trans);
    }

    print "Num bad clusters: $num_bad_clusters\n";
        
    exit($num_bad_clusters);
}

####
sub examine_cluster {
    my ($cluster, $trans_aref) = @_;

    my @trans = @$trans_aref;

    my %cluster_ids;
    foreach my $tran (@trans) {
        $tran =~ /TRINITY_(DN\d+)_/;
        $cluster_ids{$1}++;
    }

    my $num_clusters = scalar (keys %cluster_ids);
    if ($num_clusters != 1) {
        print STDERR "ERROR, got multiple clusters represented:\n"
            . "$cluster\n" . join("\n", @trans) . "\n\n";
        return(1);
    }
    else {
        return(0);
    }
}

           

