#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 TMM.fpkm.matrix trans_lengths\n\n";


my $matrix_file = $ARGV[0] or die $usage;
my $trans_lengths_file = $ARGV[1] or die $usage;

my %trans_lengths;
{
    open (my $fh, $trans_lengths_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($trans_acc, $seq_len) = split(/\s+/);
        $trans_lengths{$trans_acc} = $seq_len;
    }
    close $fh;
}




open (my $fh, $matrix_file) or die $!;
my $header = <$fh>;

my @trans;

my $sum_fpkm = 0;

while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $acc = shift @x; # gene accession
    my $max_fpkm = 0;
    my $trans_sum_fpkm = 0;
    while (@x) {
        my $fpkm = shift @x;
        
        $trans_sum_fpkm += $fpkm;
        $sum_fpkm += $fpkm;    

        if ($fpkm > $max_fpkm) {
            $max_fpkm = $fpkm;
        }
    }
    
    my $seq_len = $trans_lengths{$acc} or die "Error, no seq length for acc: $acc";
    
    push (@trans, { acc => $acc,
                    len => $seq_len,
                    sum_fpkm => $trans_sum_fpkm,
                    max_fpkm => $max_fpkm,
                });
    
}

@trans = reverse sort { $a->{sum_fpkm} <=> $b->{sum_fpkm}
                        ||
                            $a->{len} <=> $b->{len} } @trans;


## write output table
{
    open (my $ofh, ">$matrix_file.E-inputs") or die $!;
    print $ofh join("\t", "#acc", "length", "max_fpkm_over_samples", "sum_fpkm_over_samples") . "\n";
    foreach my $t (@trans) {
        print $ofh join("\t", $t->{acc}, $t->{len}, $t->{max_fpkm}, $t->{sum_fpkm}) . "\n";
    }
    close $ofh;
}


my %Estats_wanted = map { $_ => 1 } (1..100);



print "#E\tmin_fpkm\tE-N50\tnum_transcripts\n";

my $sum = 0;
my @captured;
while (@trans) {
    
    my $t = shift @trans;
    push (@captured, $t);

    $sum += $t->{sum_fpkm};

    my $pct = int($sum/$sum_fpkm * 100);
    
    if ($Estats_wanted{$pct}) {
        
        my $min_max_fpkm = &get_min_max(@captured);
        
        delete($Estats_wanted{$pct});

        my $N50 = &calc_N50(@captured);
        my $num_trans = scalar(@captured);
        
        print "E$pct\t$min_max_fpkm\t$N50\t$num_trans\n";
    }
}

# ensure that we do E100
if (%Estats_wanted) {

    my $min_max_fpkm = &get_min_max(@captured);
    my $N50 = &calc_N50(@captured);
    my $num_trans = scalar(@captured);
    
    print "E100\t$min_max_fpkm\t$N50\t$num_trans\n";
}


exit(0);


####
sub get_min_max {
    my @entries = @_;
    
    my $min = $entries[0]->{max_fpkm};
    foreach my $entry (@entries) {
        if ($entry->{max_fpkm} < $min) {
            $min = $entry->{max_fpkm};
        }
    }

    return($min);
}


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

