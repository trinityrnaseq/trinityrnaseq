#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "\n\nusage: $0 RSEM.genes.results\n\n";

my $rsem_file = $ARGV[0] or die $usage;

main: {

    my @entries;

    open (my $fh, $rsem_file) or die $!;
    my $header = <$fh>;
    my $max_fpkm = 0;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $len = $x[2];
        my $fpkm = $x[6];
        
        push (@entries, { len => $len,
                          fpkm => $fpkm,
                      });
    
        
        if ($fpkm > $max_fpkm) {
            $max_fpkm = $fpkm;
        }
    }
    close $fh;

    print join("\t", "min_fpkm", "cum_seq_len", "partial_sum_len", "N50_len", "num_entries") . "\n";
    ;
    &N50(0, \@entries);
    
    my $min_fpkm = 0.01;
    
    while ($min_fpkm < $max_fpkm) {
        
        @entries = grep { $_->{fpkm} >= $min_fpkm } @entries;
        
        &N50($min_fpkm, \@entries);
        
        $min_fpkm *= 2;
    }

}

sub N50 {
    my ($min_fpkm, $entries_aref) = @_;

    my @entries = @$entries_aref;
    @entries = reverse sort {$a->{len}<=>$b->{len}} @entries;
    
    my $num_entries = scalar(@entries);

    my $cum_seq_len = 0;
    foreach my $entry (@entries) {
        $cum_seq_len += $entry->{len};
    }
    
    my $half_cum_len = $cum_seq_len / 2;
    
    my $n50_len = "NA";
    my $partial_sum_len = 0;
    foreach my $entry (@entries) {
        my $len = $entry->{len};
        $partial_sum_len += $len;
        
        if ($partial_sum_len >= $half_cum_len) {
            $n50_len = $len;
            last;
        }
    }
    
    
    print join("\t", $min_fpkm, int($cum_seq_len+0.5), int($partial_sum_len+0.5), $n50_len, $num_entries) . "\n";

    return;
}

