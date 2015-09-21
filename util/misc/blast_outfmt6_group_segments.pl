#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use List::Util qw(min max);
use Overlap_piler;

my $usage = "\n\n\tusage: $0 blast.outfmt6 query_fasta target_fasta\n";

my $blast_file = $ARGV[0] or die $usage;
my $query_fasta = $ARGV[1] or die $usage;
my $target_fasta = $ARGV[2] or die $usage;


my %query_seq_lens = &get_seq_lengths($query_fasta);

my %target_seq_lens;
if ($query_fasta eq $target_fasta) {
    %target_seq_lens = %query_seq_lens;
}
else {
    %target_seq_lens = &get_seq_lengths($target_fasta);
}

main: {

    # header
    print join("\t", "#query_acc", "target_acc", "avg_per_id", "min_Evalue", 
               "query_match_range", "target_match_range", "pct_query_len", "pct_target_len", "max_pct_len") . "\n";


    my @hits;
    my $prev_query_target_pair = "";
    
    open (my $fh, $blast_file) or die "Error, cannot open file $blast_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $query_acc = $x[0];
        my $target_acc = $x[1];
        
        my $query_target_pair = join("$;", $query_acc, $target_acc);
        if ($query_target_pair ne $prev_query_target_pair) {
            &process_hits(@hits) if @hits;
            @hits = ();
        }
        push (@hits, [@x]);
        $prev_query_target_pair = $query_target_pair;
    }
 
    # get last one
    &process_hits(@hits);

    exit(0);
}

####
sub process_hits {
    my @hits = @_;
    
    #use Data::Dumper; print Dumper(\@hits);


    my @query_coords;
    my @target_coords;

    my $query_acc = "";
    my $target_acc = "";

    my $sum_pct_id_len = 0;
    my $sum_len = 0;
    
    my $min_evalue = 1;

    foreach my $hit (@hits) {
        my @fields = @$hit;
        
        unless ($query_acc) {
            $query_acc = $fields[0];
            $target_acc = $fields[1];
        }
        
        my ($query_lend, $query_rend) = sort {$a<=>$b} ($fields[6], $fields[7]);
        my ($target_lend, $target_rend) = sort {$a<=>$b} ($fields[8], $fields[9]);

        my $per_id = $fields[2];
        my $query_seg_len = $query_rend - $query_lend + 1;
        $sum_pct_id_len += $query_seg_len * $per_id;
        $sum_len += $query_seg_len;
        
        my $evalue = $fields[10];
        if ($evalue < $min_evalue) {
            $min_evalue = $evalue;
        }

        push (@query_coords, [$query_lend, $query_rend]);
        push (@target_coords, [$target_lend, $target_rend]);

    }
    
    my $avg_per_id = sprintf("%.2f", $sum_pct_id_len / $sum_len);
    

    my $query_match_len = 0;
    my @query_match_regions;
    {
        
        my @query_piles = &Overlap_piler::simple_coordsets_collapser(@query_coords);
        foreach my $pile (@query_piles) {
            my ($pile_lend, $pile_rend) = @$pile;
            $query_match_len += $pile_rend - $pile_lend + 1;
            push (@query_match_regions, "$pile_lend-$pile_rend");
        }
    }
    
    my $target_match_len = 0;
    my @target_match_regions;
    {
        my @target_piles = &Overlap_piler::simple_coordsets_collapser(@target_coords);
        foreach my $pile (@target_piles) {
            my ($pile_lend, $pile_rend) = @$pile;
            $target_match_len += $pile_rend - $pile_lend + 1;
            push (@target_match_regions, "$pile_lend-$pile_rend");
        }
    }
    
    my $query_len = $query_seq_lens{$query_acc} or die "Error, cannot find seq length for query: $query_acc";
    my $target_len = $target_seq_lens{$target_acc} or die "Error, cannot find seq length for target: $target_acc";

    my $pct_query_len = sprintf("%.2f", $query_match_len / $query_len * 100);
    my $pct_target_len = sprintf("%.2f", $target_match_len / $target_len * 100);
    
    
    print join("\t", $query_acc, $target_acc, $avg_per_id, $min_evalue, 
               join(",", @query_match_regions), join(",", @target_match_regions),
               $pct_query_len, $pct_target_len, max($pct_query_len, $pct_target_len) ) . "\n";

    
    return;
}


####
sub get_seq_lengths {
    my ($fasta_file) = @_;

    my %seq_lens;

    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $seq_len = length($seq_obj->get_sequence());

        # print STDERR "$acc => $seq_len\n";

        $seq_lens{$acc} = $seq_len;
    }

    return(%seq_lens);
}

