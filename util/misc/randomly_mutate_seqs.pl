#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use lib "$ENV{TRINITY_HOME}/PerlLib";
use Fasta_reader;
use List::Util qw (shuffle);
#use Math::Random;


my $usage = <<__EOUSAGE__;

##################################################################
#
# --fasta <string>              fasta filename
# 
# --subst_rate <float>          rate of single base substitutions
#
# --insert_rate <float>         rate of insertion
#
# --insert_size <int>           default: 1
#
# --delete_rate <float>         rate for deletions
#
# --delete_size <int>           default: 1
#
###################################################################


__EOUSAGE__

    ;

my $fasta_file;
my $subst_rate = 0;
my $insert_rate = 0;
my $insert_size = 1;
my $delete_rate = 0;
my $delete_size = 1;


&GetOptions (               
              'fasta=s' => \$fasta_file,
              'subst_rate=f' => \$subst_rate,
              
              'insert_rate=f' => \$insert_rate,
              'insert_size=i' => \$insert_size,

              'delete_rate=f' => \$delete_rate,
              'delete_size=i' => \$delete_size,
              
    );


unless ($fasta_file) {
    die $usage;
}
unless ($subst_rate || $insert_rate || $delete_rate) {
    die $usage;
}


main: {

    my $fasta_reader = new Fasta_reader($fasta_file);

    my %fasta_seqs = $fasta_reader->retrieve_all_seqs_hash();

    foreach my $acc (keys %fasta_seqs) {

        
        my $seq = $fasta_seqs{$acc};
        
        my @seqarray = &convert_to_seqarray($seq);
                
        my %seen;


        if ($subst_rate) {
            @seqarray = &mutate_seq(\@seqarray, \%seen, $subst_rate, 'substitution');
        }
        if ($insert_rate) {
            @seqarray = &mutate_via_insertion(\@seqarray, \%seen, $insert_rate, 'insertion');
        }
        if ($delete_rate) {
            @seqarray = &mutate_via_deletion(\@seqarray, \%seen, $delete_rate, 'deletion');
        }

        my $mutated_seq = join("", @seqarray);
        print ">$acc mutated\n$mutated_seq\n";
        
    }
    

    
    exit(0);
    
}

####
sub convert_to_seqarray {
    my ($seq) = @_;

    my @chars = split(//, uc $seq);
    
    return(@chars);
}

####
sub mutate_seq {
    my ($seqarray_aref, $seen_href, $mut_rate, $mut_type) = @_;

    my $seqlen = $#$seqarray_aref + 1;

    my $num_mutations = int($mut_rate * $seqlen + 0.5);
    
    for (1..$num_mutations) {
        
        my $pos = -1;
        do {
            $pos = int(rand($seqlen));

        } while ($seen_href->{$pos});
        
        $seen_href->{$pos} = 1;
        
        ## Substitutions
        if ($mut_type eq 'substitution') {
            
            my $char = $seqarray_aref->[$pos];
                        
            my @mut_chars = grep { $_ ne $char } qw(G A T C);
        
            my $mut_base = $mut_chars[ int(rand(length(@mut_chars))) ];

            $seqarray_aref->[$pos] = lc $mut_base;
        }

        ## Deletions
        elsif ($mut_type eq 'deletion') {
            for (my $i = $pos; $i < $pos + $delete_size && $i < $seqlen; $i++) {
                
                $seqarray_aref->[$i] = "";
                $seen_href->[$i] = 1;
            }
        }

        ## Insertions
        elsif ($mut_type eq 'insertion') {
            
            ## create insertion
            my $insertion_seq = "";
            my @bases = qw(G A T C);
            for (my $i = 0; $i < $insert_size; $i++) {
                my $base = $bases[int(rand(4))];
                $insertion_seq .= $base;
            }
            $seqarray_aref->[$pos] .= lc $insertion_seq;
        }
        else {
           confess "Error, do not understand mutation type: $mut_type ";
        }
        
    }
    

    return;
}


