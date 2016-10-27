#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 LR_blastn.outfmt6 iworm.fasta [glue_setting=2]\n\n";

my $blast_outfile = $ARGV[0] or die $usage;
my $iworm_fasta = $ARGV[1];
my $glue = $ARGV[2] || 2;

main: {

    my %LR_to_iworm_hits;
    my %iworm_to_LR_hits;

    {
        open (my $fh, $blast_outfile) or die "Error, cannot open file $blast_outfile";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my ($iworm_acc, $LR_acc, $bitscore) = ($x[0], $x[1], $x[11]);
            
            my $score = $iworm_to_LR_hits{$iworm_acc}->{$LR_acc};
            if ( (! defined $score) || $score < $bitscore) {
                $LR_to_iworm_hits{$LR_acc}->{$iworm_acc} = $bitscore;
                $iworm_to_LR_hits{$iworm_acc}->{$LR_acc} = $bitscore;
            }
        }
        close $fh;
    }

    my %LR_to_ordered_iworms;
    {
        foreach my $LR (keys %LR_to_iworm_hits) {

            my $LR_hits_href = $LR_to_iworm_hits{$LR};
            
            my @iworm_accs = reverse sort {$LR_hits_href->{$a} <=> $LR_hits_href->{$b}} keys %$LR_hits_href;

            $LR_to_ordered_iworms{$LR} = \@iworm_accs;
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

    ## Assign each iworm contig to the best-matching iworm contig linked to its best matching LR sequence.

    foreach my $iworm_acc (keys %iworm_to_LR_hits) {
        my $LR_hits_href = $iworm_to_LR_hits{$iworm_acc};
        
        my @LR_accs = reverse sort {$LR_hits_href->{$a} <=> $LR_hits_href->{$b} } keys %$LR_hits_href;
        
        my $top_LR = shift @LR_accs;
        
        my $ordered_iworms_aref = $LR_to_ordered_iworms{$top_LR};
        
        my $top_iworm_acc = $ordered_iworms_aref->[0];
        if ($iworm_acc eq $top_iworm_acc) {
            if ($#$ordered_iworms_aref > 0) {
                $top_iworm_acc = $ordered_iworms_aref->[1];
            }
            else {
                next; # no other hits
            }
        }
        
        # output record
        my $iworm_acc_A = $iworm_acc;
        my $iworm_index_A = $iworm_to_fasta_index{$iworm_acc_A};
        
        my $iworm_acc_B = $top_iworm_acc;
        my $iworm_index_B = $iworm_to_fasta_index{$iworm_acc_B};

        print join("\t", 
                   $iworm_acc_A, $iworm_index_A, 
                   $iworm_acc_B, $iworm_index_B,
                   $glue) . "\n";
    }
 

    exit(0);
}

