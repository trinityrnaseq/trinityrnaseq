#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 gmap.gff3 transcripts.fasta\n\n";

my $gmap_gff3 = $ARGV[0] or die $usage;
my $trans_fa = $ARGV[1] or die $usage;


main: {

    my %trans_coords_matched;

    open (my $fh, $gmap_gff3) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $info = $x[8];
        
        $info =~ /Target=(\S+) (\d+) (\d+)/ or die "Error, cannot parse transcript and coordinates from $info";

        my $acc = $1;
        my $lend = $2;
        my $rend = $3;
        my $per_id = $x[5];

        push (@{$trans_coords_matched{$acc}}, [$lend, $rend, $per_id]);
        
    }
    close $fh;


    my $fasta_reader = new Fasta_reader($trans_fa);
    my %trans_seqs = $fasta_reader->retrieve_all_seqs_hash();
    
    foreach my $trans_acc (keys %trans_seqs) {

        my $trans_length = length($trans_seqs{$trans_acc}) or die "Error, no trans seq for $trans_acc";

        unless (exists $trans_coords_matched{$trans_acc}) {
            print join("\t", $trans_acc, 0, $trans_length, 0, 0) . "\n";
            next;
        }
        
        ## compute length and percent identity
        my @coords = @{$trans_coords_matched{$trans_acc}};

        my $match_len = 0;
        my $sum_per_id_len = 0;
        foreach my $coordset (@coords) {
            my ($lend, $rend, $per_id) = @$coordset;
            my $len = $rend - $lend + 1;
            $match_len += $len;
            $sum_per_id_len += $len * $per_id;
        }
        my $avg_per_id = $sum_per_id_len / $match_len;
        

        my $per_length = sprintf("%.2f", $match_len / $trans_length * 100);

        print "$trans_acc\t$match_len\t$trans_length\t$per_length\t$avg_per_id\n";
    }


    exit(0);
}



