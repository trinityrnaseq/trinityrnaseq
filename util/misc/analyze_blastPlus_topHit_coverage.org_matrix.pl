#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 (count|percent) A.org_rep B.org_rep ...\n\n";

my $type = shift @ARGV;

my @org_rep_files = @ARGV or die $usage;

unless ($type eq 'count' || $type eq 'percent') { die $usage; }


main: {

    my %matrix;
    my %species_to_total_percentage;
    
    foreach my $org_rep_file (@org_rep_files) {
        open (my $fh, $org_rep_file) or die $!;
        while (<$fh>) {
            chomp;
            my ($species, $count, $percent) = split(/\t/);
            $percent =~ s/\%//;

            my $value = $percent;
            if ($type eq 'count') {
                $value = $count;
            }

            $matrix{$species}->{$org_rep_file} = $value;
            
            $species_to_total_percentage{$species} += $value;
        }
        close $fh;
    }


    my @species = reverse sort {$species_to_total_percentage{$a} <=> $species_to_total_percentage{$b}} keys %species_to_total_percentage;

    print "#\t" . join("\t", @org_rep_files) . "\n";

    foreach my $specie (@species) {
        print "$specie";
        foreach my $org (@org_rep_files) {

            my $count = $matrix{$specie}->{$org} || 0;
            print "\t$count";
        }
        print "\n";
    }
    

    exit(0);
}

