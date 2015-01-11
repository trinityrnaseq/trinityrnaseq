#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use List::Util qw(min max);
use Data::Dumper;

my $usage = "usage: $0 transcriptome.fasta polymorphRatePercentage\n\n";

my $transcriptome = $ARGV[0] or die $usage;
my $poly_rate = $ARGV[1] or die $usage;

if ($poly_rate < 1) {
    print STDERR "\n\n** WARNING: polymorphRatePercentage expects a percentage, and your input value is quite small: $poly_rate ** \n\n";
}

my %mutations;

main: {

    my $fasta_reader = new Fasta_reader($transcriptome);

    my $counter = 0;

    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
     
        my $num_snps = int($poly_rate/100 * length($sequence) + 0.5);
        
        if ($num_snps < 1) { next; } # exclude since not helpful here.
   
        &print_fasta("aleA;$acc", $sequence);

        $sequence = &mutate($sequence, $num_snps);

        &print_fasta("aleB;$acc", $sequence);
        
        $counter++;
        print STDERR "\r[$counter]    ";

    }
    print STDERR "\n\n";

    #print STDERR Dumper(\%mutations);
    
    exit(0);
}

        
####
sub print_fasta {
    my ($acc, $sequence) = @_;

    $sequence =~ s/(\S{60})/$1\n/g;

    chomp $sequence;
    
    print ">$acc\n$sequence\n";
    
    return;
}


####
sub mutate {
    my ($sequence, $num_snps) = @_;
    
    my %seen;
    
    my @chars = qw(G A T C);

    
    my @seq = split(//, uc $sequence);

    for (1..$num_snps) {
        
        my $pos;  # select unique sites (sampling w/o replacement)
    
        do {
            $pos = int(rand(length($sequence)));
            
        } while ($seen{$pos});
        
        $seen{$pos} = 1;
        
        my $nuc = uc $seq[$pos];
        
        my @others = grep { $_ ne $nuc } @chars;

        my $substitution = $others[ int(rand(scalar(@others))) ];

        #print STDERR "$nuc -> $substitution\n";

        #$mutations{"$nuc,$substitution"}++;
        
        $seq[$pos] = lc $substitution;
        
    }

    $sequence = join("", @seq);

    return($sequence);
}


            
    
