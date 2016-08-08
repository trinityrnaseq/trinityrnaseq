#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use SAM_reader;
use SAM_entry;
use Nuc_translator;

use Carp;
use Data::Dumper;


my $DEBUG = 0;

my $usage = "usage: $0 file.sam out_prefix [STRICT]\n\n";


my $sam_file = $ARGV[0] or die $usage;
my $out_prefix = $ARGV[1] or die $usage;
my $STRICT_FLAG = $ARGV[2] || 0;


my $left_fq_filename = "$out_prefix.left.fq";
my $right_fq_filename = "$out_prefix.right.fq";

open(my $left_ofh, ">$left_fq_filename") or die "Error, cannot write to $left_fq_filename";
open(my $right_ofh, ">$right_fq_filename") or die "Error, cannot write to $right_fq_filename";

main: {

    my $sam_reader = new SAM_reader($sam_file);
        
    my %fhs;

    my $prev_read_name = "";
    
    my @left_entries;
    my @right_entries;

    while (my $sam_entry = $sam_reader->get_next()) {
        
        my $core_read_name = $sam_entry->get_read_name();

        print STDERR "processing $core_read_name\n" if $DEBUG;
        
        if (! $sam_entry->is_paired()) {
            confess "ERROR, only paired reads should exist in bam file.  Encountered unpaired read: " . Dumper($sam_entry);
        }
        
        if ($prev_read_name && $core_read_name ne $prev_read_name) {
            
            print STDERR "\t-printing record for $core_read_name\n" if $DEBUG;
            
            &process_entry($prev_read_name, \@left_entries, \@right_entries);
            @left_entries = ();
            @right_entries = ();

            unless ($core_read_name gt $prev_read_name) {
                #confess "Error, it appears the sam file is not sorted by read name, as $core_read_name ! > $prev_read_name";
                # no longer die here: different tools sort in different ways, where samtools is using some mixed string and integer sorting that's not lexicographical.
            }
        }

        $prev_read_name = $core_read_name;
        if ($sam_entry->is_first_in_pair()) {
            push (@left_entries, $sam_entry);
        }
        elsif ($sam_entry->is_second_in_pair()) {
            push (@right_entries, $sam_entry);
        }
        else {
            confess "Error, cannot determine sam entry as first or second in pair: " . Dumper($sam_entry);
        }
    }
    
    ## get last one
    &process_entry($prev_read_name, \@left_entries, \@right_entries);


    exit(0);
}


####
sub process_entry {
    my ($core_read_name, $left_entries_aref, $right_entries_aref) = @_;

    unless (@$left_entries_aref) {
        print STDERR "WARNING: $core_read_name is missing first-read entry in sam file.  skipping...\n";
        if ($STRICT_FLAG) {
            confess "ERROR: $core_read_name is missing first-read entry in sam file, STRICT mode enabled";
        }
        return;
    }
    unless (@$right_entries_aref) {
        print STDERR "WARNING: $core_read_name is missing second-read entry in sam file. skipping...\n";
        if ($STRICT_FLAG) {
            confess "ERROR: $core_read_name is missing second-read entry in sam file. STRICT mode enabled";
        }
        
        return;
    }
    
    my $left_sam_entry = $left_entries_aref->[0];
    &print_fastq_record($left_ofh, $left_sam_entry);
    
    my $right_sam_entry = $right_entries_aref->[0];
    &print_fastq_record($right_ofh, $right_sam_entry);
 
    return;
}

####
sub print_fastq_record {
    my ($ofh, $sam_entry) = @_;

    
    my $read_name = $sam_entry->reconstruct_full_read_name();
    
    my $seq = $sam_entry->get_sequence();
    my $quals = $sam_entry->get_quality_scores();
    
    my $strand = $sam_entry->get_query_strand();
    if ($strand eq '-') {
        $seq = &reverse_complement($seq);
        $quals = join("", reverse(split(//, $quals)));
    }
    
    print $ofh join("\n", "\@$read_name", $seq, "+", $quals) . "\n";

    return;
}

