#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 RSEM.(isoforms|genes).results (counts|FPKM)\n\n";

my $rsem_file = $ARGV[0] or die $usage;
my $quant_type = $ARGV[1] or die $usage;

unless ($quant_type eq "counts" || $quant_type eq "FPKM") { 
    die "Error, must specify 'counts' or 'FPKM' as parameter.";
}

my $quant_column = ($quant_type eq "counts") ? 4 : 6;

my @data;

my $sum_quant = 0;

open (my $fh, $rsem_file) or die $!;
my $header = <$fh>;
while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $length = $x[2];
    my $quant = $x[$quant_column];
    push (@data, { quant => $quant,
                   length => $length,
               }
          );

    $sum_quant += $quant;
    
}
close $fh;

## sort by quant, then by length desc
@data = sort {$b->{quant}<=>$a->{quant}
              ||
                  $b->{length}<=>$a->{length}} @data;

my $sum_length = 0;
my $num_contigs = 0;
my $increment_quant = 0;
my %seen_pct;



print join("\t", "N($quant_type)", "num_contigs", "sum_contig_lengths") . "\n";
foreach my $entry (@data) {
    my $quant = $entry->{quant};
    my $length = $entry->{length};
    
    $increment_quant += $quant;
    
    my $pct_quant = sprintf("%.2f", $increment_quant/$sum_quant*100);
    
    $num_contigs++;
    $sum_length += $length;

    my $N_pct = int($pct_quant);
    unless ($seen_pct{$N_pct}) {
    
        print join("\t", $N_pct, $num_contigs, $sum_length) . "\n";
        $seen_pct{$N_pct} = 1;
    }
}

exit(0);


