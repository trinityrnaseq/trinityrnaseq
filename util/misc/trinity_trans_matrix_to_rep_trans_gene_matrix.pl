#!/usr/bin/env perl

use strict;
use warnings;

## sums up reads per gene, but uses the single isoform with the highest sum counts as the representative entry in the reported 'gene' matrix.

my $usage = "usage: $0 trinity_trans_matrix >  trinity_rep_trans_matrix\n\n";

my $trans_matrix_file = $ARGV[0] or die $usage;

main: {

    open(my $fh, $trans_matrix_file) or die $!;
    my $header = <$fh>;

    print $header;

    my %data;
    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $acc = $x[0];

        #c1000023_g1_i3

        my $gene_id;
        if ($acc =~ /^([\w_]+_g\d+)_i\d+/) {
            $gene_id = $1;
        }
        unless ($gene_id) {
            die "Error, cannot extract gene identifier for $acc";
        }

        push (@{$data{$gene_id}}, $line);
    }
    close $fh;

    foreach my $acc (keys %data) {
        my @lines = @{$data{$acc}};
        if (scalar @lines > 0) {
            my $rep_line = &get_rep(@lines);
            print "$rep_line\n";
        }
        else {
            my $rep_line = shift @lines;
            print "$rep_line\n";
        }
    }

    exit(0);
}

###
sub get_rep {
    my @lines = @_;

    my @vals;
    my $best_rowsum = -1;
    my $best_acc = "";
    
    foreach my $line (@lines) {
        my @x = split(/\t/, $line);
        my $acc = shift @x;
        my $sum = 0;
        for (my $i = 0; $i <= $#x; $i++) {
            my $val = $x[$i];
            $sum += $val;
            $vals[$i] += $val;
        }
        if ($sum > $best_rowsum) {
            $best_rowsum = $sum;
            $best_acc = $acc;
            
        }
    }
    
    my $line = join("\t", $best_acc, @vals);
    return($line);
}
    
