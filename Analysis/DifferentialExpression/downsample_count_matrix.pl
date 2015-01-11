#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(shuffle);

my $usage = "usage: $0 counts.matrix num_reads_sample\n\n";

my $matrix_file = $ARGV[0] or die $usage;
my $num_reads_to_sample = $ARGV[1] or die $usage;

main: {

    print STDERR "-parsing matrix: $matrix_file\n";
    my %matrix_in = &parse_matrix($matrix_file);  # matrix{sample}->{gene} = count;

    my %matrix_adj; # matrix{gene}->{sample} = count;



    my @samples;
    
    foreach my $sample (keys %matrix_in) {

        print STDERR "-sampling reads for sample: $sample\n";
        
        my $gene_counts_href = $matrix_in{$sample};
        
        my $sum_counts = &sum(values %$gene_counts_href);
        
        if ($sum_counts < $num_reads_to_sample) {
            print STDERR "WARNING: sample $sample only has $sum_counts reads, less than $num_reads_to_sample to sample, skipping it.\n";
            next;
        }
    
        push (@samples, $sample);
    
        my @sampled_reads = &sample_from_read_counts($gene_counts_href, $num_reads_to_sample);
        foreach my $read (@sampled_reads) {
            $matrix_adj{$read}->{$sample}++;
        }
    }
    
    print STDERR "-outputting new matrix containing downsampled reads.\n";
    ## output new matrix.
    print "\t" . join("\t", @samples) . "\n";
    foreach my $gene (keys %matrix_adj) {
        print $gene;
        foreach my $sample (@samples) {
            my $count = $matrix_adj{$gene}->{$sample} || 0;
            print "\t$count";
        }
        print "\n";
    }
    
    exit(0);

}

####
sub sample_from_read_counts {
    my ($gene_counts_href, $num_reads_to_sample) = @_;

    my @reads;
    foreach my $gene (keys %$gene_counts_href) {
        my $count = int($gene_counts_href->{$gene} + 0.5);
        for (my $i = 1; $i <= $count; $i++) {
            push (@reads, $gene);
        }
    }

    if (scalar @reads < $num_reads_to_sample) {
        die "Error, only captured " . scalar(@reads) . ", less than $num_reads_to_sample.... shouldn't happen... ";
    }

    @reads = shuffle @reads;
    @reads = @reads[0..($num_reads_to_sample-1)];
    
    return(@reads);
}


####
sub sum {
    my @vals = @_;

    my $sum = 0;
    foreach my $val (@vals) {
        $sum += $val;
    }

    return($sum);
}


####
sub parse_matrix {
    my ($matrix_file) = @_;

    open (my $fh, $matrix_file) or die "Error, cannot open file $matrix_file";
    my $header = <$fh>;
    $header =~ s/^\s+//;
    chomp $header;
    my @fields = split(/\t/, $header);
    my $adj_fields_flag = 0;

    my %matrix;

    while (<$fh>) {
        chomp;
        my @vals = split(/\t/);
        
        unless ($adj_fields_flag) {
            if (scalar(@vals) == scalar(@fields) + 1) {
                unshift(@fields, ""); # add a column for the gene name
            }
            unless (scalar @vals == scalar @fields) {
                die "Error, number of column headers is inconsistent with number of data fields.";
            }
            $adj_fields_flag = 1;
        }
        
        my $gene = $vals[0];
        for (my $i = 1; $i <= $#fields; $i++) {
            my $field = $fields[$i];
            my $val = $vals[$i];
            $matrix{$field}->{$gene} = $val;
        }
    }
    close $fh;

    return(%matrix);
}
