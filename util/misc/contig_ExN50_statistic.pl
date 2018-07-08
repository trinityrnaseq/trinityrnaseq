#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;
use File::Basename;


my $usage = "usage: $0 EXPR.matrix Trinity.fasta\n\n";


my $matrix_file = $ARGV[0] or die $usage;
my $fasta_file = $ARGV[1] or die $usage;

unless (-s $matrix_file) {
    die "Error, cannot locate matrix file: $matrix_file";
}
unless (-s $fasta_file) {
    die "Error, cannot locate fasta file: $fasta_file";
}

my %trans_lengths;
{
    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        
        my $seq_len = length($sequence);

        $trans_lengths{$acc} = $seq_len;
    }
}

open (my $fh, $matrix_file) or die $!;
my $header = <$fh>;

my @trans;

my $sum_expr = 0;

while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $acc = shift @x; # gene accession
    my $max_expr = 0;
    my $trans_sum_expr = 0;
    while (@x) {
        my $expr = shift @x;
        
        $trans_sum_expr += $expr;
        $sum_expr += $expr;    

        if ($expr > $max_expr) {
            $max_expr = $expr;
        }
    }
    
    my $seq_len = $trans_lengths{$acc} or die "Error, no seq length for acc: $acc";
    
    push (@trans, { acc => $acc,
                    len => $seq_len,
                    sum_expr => $trans_sum_expr,
                    max_expr => $max_expr,
                });
    
}

@trans = reverse sort { $a->{sum_expr} <=> $b->{sum_expr}
                        ||
                            $a->{len} <=> $b->{len} } @trans;


## write output table


my $E_file = basename($matrix_file) . ".E-inputs";
open (my $ofh, ">$E_file") or die $!;
print $ofh join("\t", "#Ex", "acc", "length", "max_expr_over_samples", "sum_expr_over_samples") . "\n";

print "Ex\tExN50\tnum_transcripts\n";

my $prev_pct = 0;
my $sum = 0;
my @captured;
while (@trans) {
    
    my $t = shift @trans;

    $sum += $t->{sum_expr};

    my $pct = int($sum/$sum_expr * 100);
    
    print $ofh join("\t", $pct, $t->{acc}, $t->{len}, sprintf("%.1f", $t->{max_expr}), sprintf("%.1f", $t->{sum_expr})) . "\n";
    
    if ($prev_pct > 0 && $pct > $prev_pct) {
        
                        
        my $N50 = &calc_N50(@captured);
        my $num_trans = scalar(@captured);
        
        print "$prev_pct\t$N50\t$num_trans\n";
    }
    
    $prev_pct = $pct;
    
    push (@captured, $t);
}

# do last one

my $N50 = &calc_N50(@captured);
my $num_trans = scalar(@captured);
print "100\t$N50\t$num_trans\n";


exit(0);


####
sub calc_N50 {
    my @entries = @_;
    
    @entries = reverse sort {$a->{len}<=>$b->{len}} @entries;
    
    my $sum_len = 0;
    foreach my $entry (@entries) {
        $sum_len += $entry->{len};
    }

    my $loc_sum = 0;
    foreach my $entry (@entries) {
        $loc_sum += $entry->{len};
        if ($loc_sum / $sum_len * 100 >= 50) {
            return($entry->{len});
        }
    }

    return(-1); # error
}

