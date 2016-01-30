#!/usr/bin/env perl 

use strict;
use warnings;
use POSIX qw (ceil);
use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 Trinity.fasta [length_bin_size=100] [out_prefix='dist']\n\n";

my $trinity_fasta = $ARGV[0] or die $usage;
my $length_bin_size = $ARGV[1] || 100;
my $out_prefix = $ARGV[2] || "dist";


main: {


    my %comp_counter;
    my %length_bin_counter;

    my $fasta_reader = new Fasta_reader($trinity_fasta);
    
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        
        my $sequence = $seq_obj->get_sequence();
        my $sequence_length = length($sequence);

        my $bin = ceil($sequence_length/$length_bin_size);
        
        $length_bin_counter{$bin}++;
        
        my $comp_name;

        if ($acc =~ /^(.*comp\d+_c\d+)/) {
            $comp_name = $1;
        }
        elsif ($acc =~ /^(.*c\d+_g\d+)/) {
            $comp_name = $1;
        }
        else {
            die "Error, cannot get component/gene id from acc: $acc";
        }
        
        $comp_counter{$comp_name}++;
    }
    
    ## get component dist
    my %comp_size_counter;
    foreach my $count (values %comp_counter) {
        $comp_size_counter{$count}++;
    }

    ## write component dist
    open (my $ofh, ">$out_prefix.comp_sizes.txt") or die $!;
    print $ofh join("\t", "transPerComp", "num_components") . "\n";
    foreach my $comp_size (sort {$a<=>$b} keys %comp_size_counter) {
        
        my $comp_size_count = $comp_size_counter{$comp_size};
        print $ofh join("\t", $comp_size, $comp_size_count) . "\n";
    }
    close $ofh;

    ## write length distribution
    open ($ofh, ">$out_prefix.trans_lengths.txt") or die $!;
    print $ofh join("\t", "#length_bin", "count_trans") . "\n";
    foreach my $length_bin (sort {$a<=>$b} keys %length_bin_counter) {
        
        my $length = $length_bin * $length_bin_size;
        
        my $count = $length_bin_counter{$length_bin};

        print $ofh join("\t", $length, $count) . "\n";
    }
    close $ofh;

    
    print STDERR "\n\nDone.  See files $out_prefix.comp_sizes.txt and $out_prefix.trans_lengths.txt\n\n";
    
    
    exit(0);
}
        
   
        
    
    
