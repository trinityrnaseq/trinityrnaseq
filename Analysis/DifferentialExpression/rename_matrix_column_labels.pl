#!/usr/bin/env perl

use strict;
use warnings;

my $usage = <<__EOUSAGE__;

###############################################
# 
#  Usage: $0 matrix.txt  samples_with_relabeling.txt
#
#
#  where samples_with_relabeling.txt has the format:
#
#   sample_name <tab> current_column_name <tab> new_column_name
#   ....
#
#   Only those entries with new names listed will be updated, so leave new_column_name blank for those that should remain unchanged.
#
#
#################################################


__EOUSAGE__

    ;

my $matrix_file = $ARGV[0] or die $usage;
my $samples_with_relabeling = $ARGV[1] or die $usage;

main: {

    my %new_column_names;
    {
        open (my $fh, $samples_with_relabeling) or die $!;
        while (<$fh>) {
            chomp;
            my ($sample, $old_name, $new_name) = split(/\t/);
            if ($new_name) {
                $new_column_names{$old_name} = $new_name;
            }
            
        }
        close $fh;
    }

    open (my $fh, $matrix_file) or die $!;
    my $header = <$fh>;
    chomp $header;
    my @fields = split(/\t/, $header);
    
    foreach my $field (@fields) {
        if (my $new_name = $new_column_names{$field}) {
            $field = $new_name;
        }
    }
    print join("\t", @fields) . "\n";
    while (<$fh>) {
        print $_;
    }
    

    exit(0);
}
