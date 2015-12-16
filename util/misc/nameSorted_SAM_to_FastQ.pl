#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use SAM_reader;
use SAM_entry;
use Nuc_translator;

my $usage = "usage: $0 file.sam out_prefix\n\n";


my $sam_file = $ARGV[0] or die $usage;
my $out_prefix = $ARGV[1] or die $usage;



main: {

    my $sam_reader = new SAM_reader($sam_file);
        
    my %fhs;
 
    while (my $sam_entry = $sam_reader->get_next()) {
        
        if ($sam_entry->is_query_unmapped()) {
            next;
        }

        my $outfh = &get_output_file($sam_entry, \%fhs);

        my $seq = $sam_entry->get_sequence();
        my $quals = $sam_entry->get_quality_scores();

        my $strand = $sam_entry->get_query_strand();
        if ($strand eq '-') {
            $seq = &reverse_complement($seq);
            $quals = join("", reverse(split(//, $quals)));
        }
        
        my $read_name = $sam_entry->reconstruct_full_read_name();

        print $outfh join("\n", "\@$read_name", $seq, "+", $quals) . "\n";
    }

    exit(0);
}
        
####
sub get_output_file {
    my ($sam_entry, $fhs_href) = @_;

    my $filename = "$out_prefix.single.fq";

    if ($sam_entry->is_paired()) {
        
        if ($sam_entry->is_first_in_pair()) {

            $filename = "$out_prefix.left.fq";
        }
        elsif ($sam_entry->is_second_in_pair()) {

            $filename = "$out_prefix.right.fq";
        }
        else {
            die "Error, read is paired but neither first or second in pair!";
        }
    }

    my $fh = $fhs_href->{$filename};

    unless ($fh) {
        
        open ($fh, ">$filename") or die "Error, cannot write to $filename";
        $fhs_href->{$filename} = $fh;
    }


    return($fh);
}


