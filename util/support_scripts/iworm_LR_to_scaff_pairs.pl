#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

my $usage = "\n\n\tusage: $0 LR_blastn.outfmt6 iworm.fasta min_containment=98 KMER_LENGTH=25\n\n";

my $blast_outfile = $ARGV[0] or die $usage;
my $iworm_fasta = $ARGV[1];
my $MIN_CONTAINMENT = $ARGV[2] || 98;
my $KMER_LENGTH = $ARGV[3] || 25;

main: {

    my %LR_to_iworm_hits;

    {
        open (my $fh, $blast_outfile) or die "Error, cannot open file $blast_outfile";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my ($iworm_acc, $LR_acc, 
                $iworm_end5, $iworm_end3,
                $LR_end5, $LR_end3,
                $bitscore,
                $max_either_containment,
                ) = ($x[0], $x[1], 
                     $x[6], $x[7], 
                     $x[8], $x[9],
                     $x[11],
                     $x[16],
                    );
            
            unless ($max_either_containment >= $MIN_CONTAINMENT) { next; }
            

            my ($LR_lend, $LR_rend) = sort {$a<=>$b} ($LR_end5, $LR_end3);

            push (@{$LR_to_iworm_hits{$LR_acc}}, { iworm_acc => $iworm_acc,
                                                   lend => $LR_lend,
                                                   rend => $LR_rend,
                  });
            
            
            
        }
        close $fh;
    }

    
    my %iworm_pairs;
    
    foreach my $LR (keys %LR_to_iworm_hits) {

        my @iworm_hits = sort {$a->{lend}<=>$b->{lend}
                               ||
                                   $a->{rend} <=> $b->{rend}} @{$LR_to_iworm_hits{$LR}};

        
        for (my $i = 0; $i < $#iworm_hits; $i++) {

            my $iworm_i = $iworm_hits[$i];

            for (my $j = $i + 1; $j <= $#iworm_hits; $j++) {

                my $iworm_j = $iworm_hits[$j];

                if ($iworm_j->{lend} > $iworm_i->{rend}) {
                    last; 
                }

                if (&overlap_by_Kminus1($iworm_i, $iworm_j, $KMER_LENGTH)) {
                    
                    my $pair_token = join("$;", sort ($iworm_i->{iworm_acc}, $iworm_j->{iworm_acc}));
                    $iworm_pairs{$pair_token}++;
                    
                }
            }
        }

    }
    

    ## get fasta index values for the iworm contigs:
    my %iworm_to_fasta_index;
    {
        my $index = 0;
        open (my $fh, $iworm_fasta) or die "Error, cannot open file $iworm_fasta";
        while (<$fh>) {
            chomp;
            if (/^>(\S+)/) {
                my $acc = $1;
                $iworm_to_fasta_index{$acc} = $index;
                $index++;
            }
        }
        close $fh;
    }

    
    foreach my $pair (keys %iworm_pairs) {

        my ($iworm_acc_A, $iworm_acc_B) = split(/$;/, $pair);
        
        # a179572;19
        
        # assign glue as the lower of the coverage values:

        my ($coreA, $covA) = split(/;/, $iworm_acc_A);
        my ($coreB, $covB) = split(/;/, $iworm_acc_B);
        
        my $glue = min($covA, $covB);
        
        # output record
        my $iworm_index_A = $iworm_to_fasta_index{$iworm_acc_A};
        
        my $iworm_index_B = $iworm_to_fasta_index{$iworm_acc_B};

                
        print join("\t", 
                   $iworm_acc_A, $iworm_index_A, 
                   $iworm_acc_B, $iworm_index_B,
                   $glue) . "\n";
    }



    exit(0);
    
}


####
sub overlap_by_Kminus1 {
    my ($iworm_i, $iworm_j, $K) = @_;

    ## require:
    ##
    ##   -----------------   iworm_i
    ##             -------------------   iworm_j
    ##
    ## ----------------------------------  LR
    ##
    ##  overlap by at least K-1
    ##  staggered overlaps along LR
    ##  iworm_i and iworm_j alignment lengths >= 2 * (K-1)

    
    my $iworm_i_align_len = $iworm_i->{rend} - $iworm_i->{lend} + 1;
    my $iworm_j_align_len = $iworm_j->{rend} - $iworm_j->{lend} + 1;

    my $min_align_len = 2 * ($K-1);
    
    if ($iworm_i_align_len < $min_align_len 
        ||
        $iworm_j_align_len < $min_align_len) { 
       
        return(0);
    }

    if ($iworm_i->{rend} > $iworm_j->{rend}) {
        return(0);
    }
    
    
    my $overlap_len = $iworm_i->{rend} - $iworm_j->{lend} + 1;
    
    if ($overlap_len >= $K - 1) {
        
        return(1);
    }
    else {
        return(0);
    }
    
}

