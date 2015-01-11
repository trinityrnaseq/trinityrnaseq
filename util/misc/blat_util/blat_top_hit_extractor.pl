#!/usr/bin/env perl

use strict;

my $SEE = 0;

###############
# blat format: # Q=cDNA T=genomic
################

#  0: match
#  1: mis-match
#  2: rep. match
#  3: N's
#  4: Q gap count
#  5: Q gap bases
#  6: T gap count
#  7: T gap bases
#  8: strand
#  9: Q name
# 10: Q size
# 11: Q start
# 12: Q end
# 13: T name
# 14: T size
# 15: T start
# 16: T end
# 17: block count
# 18: block Sizes
# 19: Q starts
# 20: T starts
# 21: Q seqs (pslx format)
# 22: T seqs (pslx format)

#############################################################
## This script filters blat output and extracts only       ##
## the top scoring alignment chain for each accession.     ##
#############################################################

my $line_num = 0;
my %data;

my $filename = $ARGV[0] or die "\n\nusage: $0 outputfile.psl [num_top_hits=1]\n\n";
my $num_top_hits = $ARGV[1] || 1;

## First Pass, assign scores to entries associated with accessions.
open (FILE, "$filename");
my $line_num = 0;
while (<FILE>) {
    $line_num++;
    my @x = split (/\t/);
    unless ($x[0] =~ /^\d/) {next;}
    my ($matches, $mismatches, $q_gap_num, $t_gap_num, $q_insert, $t_insert) = ($x[0], $x[1], $x[4], $x[6], $x[5], abs($x[7]));
    
    my $score = &calculate_score($matches, $mismatches, $q_gap_num, $t_gap_num, $q_insert, $t_insert);
    
    my $accession = $x[9];
    my $score_struct = {score=>$score,
                        line_num=>$line_num};
    push (@{$data{$accession}}, $score_struct);
}

close FILE;

## Identify each highest scoring match
my %line_nums_to_print;

foreach my $accession (keys %data) {
    my @hits = @{$data{$accession}};
    @hits = reverse sort {$a->{score}<=>$b->{score}} @hits; #sort in reverse order of score.
    
    if ($SEE) { #verify the top hit is chosen.
        foreach my $hit (@hits) {
            print $hit->{line_num} . ":" . $hit->{score} . " ";
        }
        
        print "\n chose ";
    }
    
    for (my $i = 0; $i < $num_top_hits && $i <= $#hits; $i++) {
        
        my $top_hit = $hits[$i];
        print $top_hit->{line_num} . ":" . $top_hit->{score} . "\n" if $SEE;
        my $line_num = $top_hit->{line_num};
        $line_nums_to_print{$line_num} = 1;
    }
}

open (FILE, "$filename");
$line_num = 0;
while (<FILE>) {
    $line_num++;
    if ($line_nums_to_print{$line_num}) {
        print;
    }
}
close FILE;

exit(0);


####
sub calculate_score { 
    my ($match, $mismatch, $q_gap_num, $t_gap_num, $q_insert, $t_insert) = @_;
    
    ## JKent's code in pslFilter:
    # score = (psl->match + psl->repMatch)*reward - psl->misMatch*cost 
    #    - (psl->qNumInsert + psl->tNumInsert + 1) * gapOpenCost
    #	- log(psl->qBaseInsert + psl->tBaseInsert + 1) * gapSizeLogMod;
    
    #use jkent's default score parameters:
    my $reward = 1;
    my $cost = 1;
    my $gapOpenCost = 4;
    my $gapSizeLogMod = 1;
    
    return ( ($match * $reward) - ($mismatch * $cost) - 
             ( ($q_gap_num + $t_gap_num) * $gapOpenCost) - 
             ( log ($q_insert + $t_insert + 1) * $gapSizeLogMod) );
    
} 
