#!/usr/bin/env perl

use strict;
use warnings;


################
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

## All sequences start at 0 here; array-based.

my $JOIN_GAP = 9; #join alignment segments if within this gap along the genomic sequence.

my $chain_number = 0;

while (<STDIN>) {
    unless (/\w/) { next; }
    chomp; 
    my @x = split (/\t/);
    unless ($x[0] =~ /^\d/) {next;} #eliminate headers if present.
    my @alignment_segments;
    $chain_number++;
    my $strand = $x[8];
    my $cdna_name = $x[9];
    my $genomic_name = $x[13];
    my $genomic_length = $x[14];
    my $cdna_length = $x[10];
    my $cDNA_seqs = $x[21];
    my $genomic_seqs = $x[22];
    my $num_segs = $x[17];
    
    my @cdna_coords = split (/,/, $x[19]);
    my @genomic_coords = split (/,/, $x[20]);
    my @lengths = split (/,/, $x[18]);
    my @cDNA_seqs = split (/,/, $x[21]);
    my @genomic_seqs = split (/,/, $x[22]);
    
    ## going to implement score as chain_score + segment score
    ## Chain score = matches_num - mismatch_num
    my $chain_score = $x[0] - $x[1];
    
    ## report each segment match as a separate btab entry:
    my $segment_number = 0;
    for (my $i = 0; $i < $num_segs; $i++) {
        $segment_number++;
        my $length = $lengths[$i];
        my $cdna_coord = $cdna_coords[$i];
        my $genomic_coord = $genomic_coords[$i];
        my ($cdna_end5, $cdna_end3) = (++$cdna_coord, $cdna_coord + $length - 1);
        my ($genomic_end5, $genomic_end3) = (++$genomic_coord, $genomic_coord + $length -1);
        if ($strand eq "-") {
            ($genomic_end5, $genomic_end3) = ($genomic_end3, $genomic_end5);
            ($cdna_end5, $cdna_end3) = sort {$a<=>$b} ($cdna_length - $cdna_end5 + 1, $cdna_length - $cdna_end3 + 1);
        }
        my $segment_score = $chain_score + $length;
        my $per_id = -1;
        my ($gseq, $cseq);
        if ( ($gseq = $genomic_seqs[$i]) && ($cseq = $cDNA_seqs[$i])) {
            if ($gseq eq $cseq) {
                $per_id = 100;
            } else {
                ## walk thru and determine num ids
                my $num_id = 0;
                my @gseq_array = split (//, $gseq);
                my @cseq_array = split (//, $cseq);
                for (my $j = 0; $j < $length; $j++) {
                    if ($gseq_array[$j] eq $cseq_array[$j]) {
                        $num_id++;
                    }
                }
                $per_id = ($num_id/$length) * 100;
            }
        }
        
        ## Create btab line.
        my @btab;
        $btab[0] = $genomic_name;
        $btab[2] = $genomic_length;
        $btab[3] = "blat";
        $btab[5] = $cdna_name;
        $btab[6] = $genomic_end5;
        $btab[7] = $genomic_end3;
        $btab[8] = $cdna_end5;
        $btab[9] = $cdna_end3;
        
        $btab[10] = $per_id;
        $btab[12] = $segment_score;
        
        $btab[13] = $chain_number;
        $btab[14] = $segment_number;
        
        $btab[18] = $length;
        push (@alignment_segments, [@btab]);
        
    }
    
    &process_alignment_chain(\@alignment_segments, $chain_number);
}

exit(0);


## Join alignment segments if within 5 bp along the genomic sequence
sub process_alignment_chain {
    my $alignment_segments_aref = shift;
    my $chain_number = shift;
    my @segments = sort {$a->[8]<=>$b->[8]} @$alignment_segments_aref;
    if ($#segments > 0) {
        my @new_segments = ($segments[0]); # always holds the last segment analyzed.
        for (my $i=1; $i <= $#segments; $i++) {
            my $last_segment = $new_segments[$#new_segments];
            my $current_segment = $segments[$i];
            my $prev_end3 = $last_segment->[7];
            my $curr_end5 = $current_segment->[6];
            my $gap_length = abs ($prev_end3 - $curr_end5) - 1;
            
            if ($gap_length <= $JOIN_GAP) {
                ## Must join prev and current segment
                my $prev_seg_length = abs ($last_segment->[7] - $last_segment->[6]) + 1;
                my $curr_seg_length = abs ($current_segment->[7] - $current_segment->[6]) + 1;
                
                ## make prev end3 the new end3 for both genomic and cdna coordinates
                $last_segment->[7] = $current_segment->[7];
                $last_segment->[9] = $current_segment->[9];
                
                ## recalculate the percent ID
                my $prev_per_id = $last_segment->[10];
                my $curr_per_id = $current_segment->[10];
                
                my $new_per_id = ($prev_seg_length * $prev_per_id + $curr_seg_length * $curr_per_id) / 
                    ($prev_seg_length + $curr_seg_length + $gap_length);
                
                $last_segment->[10] = $new_per_id;
                $last_segment->[18] = abs($last_segment->[9] - $last_segment->[8]) + 1;
                $last_segment->[12] += $curr_seg_length;
                
            } else {
                
                #make the current segment the last segment
                push (@new_segments, $current_segment);
                
            }
        }
        @segments = @new_segments;
        
        ## renumber the segment numbers
        my $segnum = 0;
        foreach my $segment (@segments) {
            $segnum++;
            $segment->[14] = $segnum;
        }
    }
    
    ## write GFF format.
    my $orient;
    foreach my $segment (@segments) {
        
        my $genomic_end5 = $segment->[6];
        my $genomic_end3 = $segment->[7];
        
        unless ($orient) {
            if ($genomic_end5 < $genomic_end3) {
                $orient = '+';
            }
            elsif ($genomic_end3 < $genomic_end5) {
                $orient = '-';
            }
        }
        if ($orient) { 
            last;
        }
    }
    
    unless ($orient) {
        $orient = '+'; ## set a default
    }


    foreach my $segment (@segments) {
        
        my $genomic_contig = $segment->[0];
        my $cdna_name = $segment->[5];
        
        my $genomic_end5 = $segment->[6];
        my $genomic_end3 = $segment->[7];
        
        my ($genomic_lend, $genomic_rend) = sort {$a<=>$b} ($genomic_end3, $genomic_end5);
        
        my $cdna_end5 = $segment->[8];
        my $cdna_end3 = $segment->[9];
        
        my $per_id = sprintf("%.2f", $segment->[10]);
        
        my $chain_ID = "blat.proc$$.chain_" . $chain_number;
        
        print join("\t", $genomic_contig, "BLAT", "cDNA_match",
                   $genomic_lend, $genomic_rend, $per_id, $orient, ".",
                   "ID=$chain_ID;Target=$cdna_name $cdna_end5 $cdna_end3 +") . "\n";
        
        
        
        
    }

    return;
}











