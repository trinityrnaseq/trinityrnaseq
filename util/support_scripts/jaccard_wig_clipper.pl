#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");

use WigParser;

my $jaccard_file;
my $trough_win = 200;
my $min_jaccard_delta = 0.35;
my $max_jaccard_trough_val = 0.05;
my $coverage_wig;

my $usage = <<__EOUSAGE__;

##############################################################################################################################
#
# Required;
#
#  --jaccard_wig <string>             :jaccard wig file
#
# Optional:
#
#  --trough_win <int>                 :default($trough_win)
#  --min_jaccard_delta <float>        :default($min_jaccard_delta)
#  --max_jaccard_trough_val <float>   :default($max_jaccard_trough_val)
#
#  --coverage_wig <string>            :fragment coverage wig. Smallest coverage value in trough_win of putative clip is used.
#
#
##############################################################################################################################


__EOUSAGE__

    ;




my $VERBOSE = 0;

my $help_flag;

&GetOptions ( 'h' => \$help_flag,
              'jaccard_wig=s' => \$jaccard_file,
              
              'trough_win=i' => \$trough_win,
              'min_jaccard_delta=f' => \$min_jaccard_delta,
              'max_jaccard_trough_val=f' => \$max_jaccard_trough_val,

              'coverage_wig=s' => \$coverage_wig,
              
              );


unless ($jaccard_file) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, couldn't parse params: @ARGV";
}


main: {


    my $scaff_to_jaccard_parser = new WigParser($jaccard_file);
    
    my $scaff_to_frag_coverage_parser;
    
    if ($coverage_wig) {
        $scaff_to_frag_coverage_parser = new WigParser($coverage_wig);
    }
    
	foreach my $scaff (sort $scaff_to_jaccard_parser->get_contig_list()) {
		
        my @jaccard_vals = $scaff_to_jaccard_parser->get_wig_array($scaff);
        
        my $vals_aref = \@jaccard_vals; # retaining old usage after module refactoring
        

        my @trough_positions;

        for (my $i = $trough_win; $i <= $#$vals_aref; $i++) {

            my $left_win_pos = $i - $trough_win;
            my $right_win_pos = $i;
            my $center_win = int( ($left_win_pos + $right_win_pos)/2);

            my $jaccard_val_left = $vals_aref->[$left_win_pos];
            my $jaccard_val_right = $vals_aref->[$right_win_pos];
            my $jaccard_val_center = $vals_aref->[$center_win];

            
            if ($jaccard_val_center <= $max_jaccard_trough_val) {
                
                ## found a trough
                push (@trough_positions, { pos => $center_win,
                                           jaccard => $jaccard_val_center,
                                           avg_delta => (  ($jaccard_val_left - $jaccard_val_center ) + ($jaccard_val_right - $jaccard_val_center) ) / 2,
                                       } );
            }

        }
        
        if (@trough_positions) {
            @trough_positions = &group_trough_positions_within_window_select_best_clip(@trough_positions);
            
            ## restrict clips to those that are in troughs of required depth from neighboring hills
            my @clips = &require_neighboring_hills(\@trough_positions, $vals_aref);
            
            if ($coverage_wig) {
                
                my @coverage_array = $scaff_to_frag_coverage_parser->get_wig_array($scaff);

                &reposition_clips_by_min_coverage(\@clips, \@coverage_array);
            }
            
            
            print "variableStep chrom=$scaff\n";
            foreach my $clip (@clips) {
                my $pos = $clip->{pos};
                
                print join("\t", $pos, 1) . "\n";
                
            }
        }
		
	}
    
	
	exit(0);
    
    
}


####
sub group_trough_positions_within_window_select_best_clip {
    my @trough_positions = @_;

    my @trough_groups;
    
    my $trough_pos = shift @trough_positions;
    push (@trough_groups, [$trough_pos]);
    
    while (@trough_positions) {
        
        my $curr_trough_entry = shift @trough_positions;
        my $curr_trough_pos = $curr_trough_entry->{pos};


        my $prev_trough_group_aref = $trough_groups[$#trough_groups];
        my $last_pos_val = $prev_trough_group_aref->[$#$prev_trough_group_aref]->{pos};
        
        if ($curr_trough_pos - $last_pos_val <= $trough_win) {

            push (@$prev_trough_group_aref, $curr_trough_entry);
        }
        else {
            ## start a new group
            push (@trough_groups, [ $curr_trough_entry ] );
        }
    }


    ## take the best one as a clip point
    ## best defined here as the one with the least jaccard support and greatest delta value
    
    my @clips;
    
    foreach my $group (@trough_groups) {
        
        my @eles = @$group;

        @eles = sort {$a->{jaccard}<=>$b->{jaccard}
                      ||
                          $b->{avg_delta} <=> $a->{avg_delta}    ## order should be low jaccard, high avg_delta
                  } @eles;

        my $clip_ele = shift @eles;
        
        push (@clips, $clip_ele);
    }


    return(@clips);
}
        
            
####
sub require_neighboring_hills {
    my ($clips_aref, $jaccard_aref) = @_;

    my @validated_clips;


    foreach my $clip (@$clips_aref) {
        my $pos = $clip->{pos};
        my $jaccard = $clip->{jaccard};
        
        my $hill_jaccard_min = $jaccard + $min_jaccard_delta;
        

        my $left_hill_stop = $pos - $trough_win;
        if ($left_hill_stop < 1) {
            $left_hill_stop = 1;
        }
        
        my $right_hill_stop = $pos + $trough_win;
        if ($right_hill_stop > $#$jaccard_aref) {
            $right_hill_stop = $#$jaccard_aref;
        }

        if (&search_for_hill($jaccard_aref, $left_hill_stop, $pos-1, $hill_jaccard_min)
            &&
            &search_for_hill($jaccard_aref, $pos + 1, $right_hill_stop, $hill_jaccard_min) ) {

            push (@validated_clips, $clip);
        }
    }


    return(@validated_clips);
}
 
####
sub search_for_hill {
    my ($jaccard_aref, $lend, $rend, $min_jaccard_val) = @_;

    for (my $i = $lend; $i <= $rend; $i++) {
        if ($jaccard_aref->[$i] >= $min_jaccard_val) {
            return(1);
        }
    }

    return(0);
}

      
####
sub reposition_clips_by_min_coverage {
    my ($clips_aref, $coverage_aref) = @_;

    foreach my $clip (@$clips_aref) {

        my $clip_pos = $clip->{pos};
        
        my $win_left = $clip_pos - int($trough_win/2);
        if ($win_left < 1) {
            $win_left = 1;
        }
        my $win_right = $clip_pos + int($trough_win/2);
        if ($win_right > $#$coverage_aref) {
            $win_right = $#$coverage_aref;
        }


        my $min_cov_val = $coverage_aref->[$clip_pos];
        my $min_cov_pos = $clip_pos;
        
        for (my $i = $win_left; $i <= $win_right; $i++) {
            my $cov = $coverage_aref->[$i];
            if ($cov < $min_cov_val) {
                $min_cov_val = $cov;
                $min_cov_pos = $i;
            }
        }
        $clip->{pos} = $min_cov_pos;
    }

    return;
}


        
