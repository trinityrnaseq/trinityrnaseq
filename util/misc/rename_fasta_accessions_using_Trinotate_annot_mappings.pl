#!/usr/bin/env perl

use strict;
use warnings;

my $usage = <<__EOUSAGE__;

###############################################
# 
#  Usage: $0 Trinity.fasta  new_feature_id_mapping.txt
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

my $trinity_fasta_file = $ARGV[0] or die $usage;
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
    
    open (my $fh, $trinity_fasta_file) or die $!;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)/) {
            my $acc = $1;
            if (my $new_id = $new_ids{$acc}) {
                $line =~ s/^>/^>$new_id /;
                #print $line;
            }
        }
        print $line;
    }
    

    exit(0);
}
