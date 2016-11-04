#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 scaff_pairs_A.txt scaff_pairs_B.txt ...\n\n";

my @pair_files = @ARGV;
unless (scalar(@pair_files) > 1) {
    die $usage;
}

main: {

    my %top_pairs;
    foreach my $pair_file (@pair_files) {
        open (my $fh, $pair_file) or die "Error, cannot open file $pair_file";
        while (<$fh>) {
            my $line = $_;
            chomp;
            
            my ($iworm_A, $idx_A, $iworm_B, $idx_B, $count) = split(/\t/);

            my $pair_token = join("$;", sort ($iworm_A, $iworm_B));

            if ( (! exists $top_pairs{$pair_token})
                 ||
                 $top_pairs{$pair_token}->{count} < $count) {

                $top_pairs{$pair_token} = { count => $count,
                                            line => $line };

            }
        }
        close $fh;
    }


    my @top_structs = reverse sort {$a->{count} <=> $b->{count}} values %top_pairs;

    foreach my $struct (@top_structs) {

        print $struct->{line};
    
    }

    exit(0);
}

