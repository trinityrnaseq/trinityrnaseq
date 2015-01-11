#!/usr/bin/env perl

use strict;
use warnings;

my $usage = <<__EOUSAGE__;

###############################################
# 
#  Usage: $0 matrix.txt  new_feature_id_mapping.txt
#
#  The 'new_feature_id_mapping.txt' file has the format:
#
#   current_identifier <tab> new_identifier
#   ....
#
#
#   Only those entries with new names listed will be updated, the rest stay unchanged.
#
#
#################################################


__EOUSAGE__

    ;

my $matrix_file = $ARGV[0] or die $usage;
my $new_feature_mappings = $ARGV[1] or die $usage;

main: {

    my %new_ids;
    {
        open (my $fh, $new_feature_mappings) or die $!;
        while (<$fh>) {
            chomp;
            my ($old_name, $new_name) = split(/\t/);
            $new_ids{$old_name} = $new_name;
        }
        close $fh;
    }
    
    open (my $fh, $matrix_file) or die $!;
    while (<$fh>) {
        unless(/\w/) {
            print;
            next;
        }
        
        chomp;
        my @ids = split(/\t/);
        if (my $new_id = $new_ids{$ids[0]}) {
            $ids[0] = $new_id;
        }
        if (scalar(@ids) > 1 && (my $new_id = $new_ids{$ids[1]})) {
            $ids[1] = $new_id;
        }
        
        print join("\t", @ids) . "\n";
    }
    

    exit(0);
}
