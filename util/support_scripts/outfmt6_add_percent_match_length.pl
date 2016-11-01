#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use List::Util qw(min max);

my $usage = "\n\n\tusage: $0 blast.outfmt6 query_fasta target_fasta\n";

my $blast_file = $ARGV[0] or die $usage;
my $query_fasta = $ARGV[1] or die $usage;
my $target_fasta = $ARGV[2] or die $usage;

main: {

    my %query_seq_lens = &get_seq_lengths($query_fasta);

    my %target_seq_lens;
    if ($query_fasta eq $target_fasta) {
        %target_seq_lens = %query_seq_lens;
    }
    else {
        %target_seq_lens = &get_seq_lengths($target_fasta);
    }

    open (my $fh, $blast_file) or die "Error, cannot open file $blast_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $query_acc = $x[0];
        my $target_acc = $x[1];

        my $query_len = $query_seq_lens{$query_acc} or die "Error, cannot find seq length for query: $query_acc";
        my $target_len = $target_seq_lens{$target_acc} or die "Error, cannot find seq length for target: $target_acc";

        my $query_hit_len = abs($x[7]-$x[6]);
        my $db_hit_len = abs($x[9]-$x[8]);

        my $pct_query_len = sprintf("%.2f", $query_hit_len / $query_len * 100);
        my $pct_target_len = sprintf("%.2f", $db_hit_len / $target_len * 100);

        push (@x, $query_len, $pct_query_len, $target_len, $pct_target_len, max($pct_query_len, $pct_target_len));
        
        print join("\t", @x) . "\n";
    }
    
    exit(0);
}

####
sub get_seq_lengths {
    my ($fasta_file) = @_;

    my %seq_lens;

    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $seq_len = length($seq_obj->get_sequence());

        $seq_lens{$acc} = $seq_len;
    }

    return(%seq_lens);
}

