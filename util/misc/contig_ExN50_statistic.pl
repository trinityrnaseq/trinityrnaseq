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

my %gene_to_trans;

my $sum_expr = 0;

my $feature_type = "transcript";

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

    my $gene_id = $acc;
    if ($acc =~ /^(\S+)_i\d+/) {
        $gene_id = $1;
        $feature_type = "gene";
    }
    
    push (@{$gene_to_trans{$gene_id}}, { acc => $acc,
                                         len => $seq_len,
                                         sum_expr => $trans_sum_expr,
                                         max_expr => $max_expr,
          });
    
}

my @genes;

## make expression weighted gene length
foreach my $gene (keys %gene_to_trans) {
    my @trans_structs = @{$gene_to_trans{$gene}};

    my $sum_expr = 0;
    my $sum_expr_n_len = 0;
    my $max_expr = 0;
    foreach my $trans_struct (@trans_structs) {
        my $len = $trans_struct->{len};
        my $expr = $trans_struct->{sum_expr} || 1;
        $sum_expr_n_len += $len * $expr;
        $sum_expr += $expr;

        my $m_expr = $trans_struct->{max_expr};
        if ($m_expr > $max_expr) {
            $max_expr = $m_expr;
        }
        
    }
    
    my $gene_len = $sum_expr_n_len / $sum_expr;
    
    push (@genes, { acc => $gene,
                    sum_expr => $sum_expr,
                    len => $gene_len,
                    max_expr => $max_expr } );
}


@genes = reverse sort { $a->{sum_expr} <=> $b->{sum_expr}
                        ||
                            $a->{len} <=> $b->{len} } @genes;



## write output table


my $E_file = basename($matrix_file) . ".E-inputs";
open (my $ofh, ">$E_file") or die $!;
print $ofh join("\t", "#Ex", "acc", "length", "max_expr_over_samples", "sum_expr_over_samples") . "\n";

print "Ex\tExN50\tnum_${feature_type}s\n";

my $prev_pct = 0;
my $sum = 0;
my @captured;
while (@genes) {
    
    my $t = shift @genes;

    $sum += $t->{sum_expr};

    my $pct = int($sum/$sum_expr * 100);
    
    print $ofh join("\t", $pct, $t->{acc}, int($t->{len}), sprintf("%.1f", $t->{max_expr}), sprintf("%.1f", $t->{sum_expr})) . "\n";
    
    if ($prev_pct > 0 && $pct > $prev_pct) {
        
                        
        my $N50 = int(&calc_N50(@captured));
        my $num_trans = scalar(@captured);
        
        print "$prev_pct\t$N50\t$num_trans\n";
    }
    
    $prev_pct = $pct;
    
    push (@captured, $t);
}

# do last one

my $N50 = &calc_N50(@captured);
my $num_genes = scalar(@captured);
print "100\t$N50\t$num_genes\n";


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

