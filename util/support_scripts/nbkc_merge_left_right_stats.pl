#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use DelimParser;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

###############################################################################
#
# Required:
#
#  --left <string>     left.fq.stats 
#  --right <string>    right.fq.stats
#
# Optional
#
#  --sorted            flag indicating that entries are lexically sorted
#                       (this can account for differences in representation by 
#                        reads in either file)
#                       Unpaired entries are ignored.
#
################################################################################

__EOUSAGE__

    ;

my $left_stats_file;
my $right_stats_file;

my $sorted_flag = 0;

&GetOptions( 'left=s' => \$left_stats_file,
             'right=s' => \$right_stats_file,
             'sorted' => \$sorted_flag,
             );


unless ($left_stats_file && $right_stats_file) {
    die $usage;
}


main: {

    my ($left_fh, $right_fh);
    print STDERR "-opening $left_stats_file\n";
    if ($left_stats_file =~ /\.gz$/) {
        open ($left_fh, "gunzip -c $left_stats_file | ") or die $!;
    } elsif ($left_stats_file =~ /\.xz$/) {
        open(${left_fh}, "xz -cd ${left_stats_file} | ") or die $!;
    }
    else {
        open ($left_fh, $left_stats_file) or die $!;
    }

    print STDERR "-opening $right_stats_file\n";
    if ($right_stats_file =~ /\.gz$/) {
        open ($right_fh, "gunzip -c $right_stats_file | ") or die $!;
    } elsif ($right_stats_file =~ /\.xz$/) {
        open (${right_fh}, "xz -dc ${right_stats_file} | ") or die $!;
    }
    else {
        open ($right_fh, $right_stats_file) or die $!;
    }

    print STDERR "-done opening files.\n";
    
    my $left_delim_parser = new DelimParser::Reader($left_fh, "\t");
    my $right_delim_parser = new DelimParser::Reader($right_fh, "\t");

    
    # print header:
    print join("\t", "acc", 
               "left_acc", "left_median_cov", "left_mean_cov", "left_stdev",
               "right_acc", "right_median_cov", "right_mean_cov", "right_stdev",
               "median_cov", "mean_cov", "stdev"
        ) . "\n";
    
    my $left_row = $left_delim_parser->get_row();
    my $right_row = $right_delim_parser->get_row();

    my $num_left = ($left_row) ? 1 : 0;
    my $num_right = ($right_row) ? 1 : 0;
    my $num_paired = 0;
    my ($example_left_acc, $example_right_acc) = (($left_row) ? $left_row->{acc} : "", ($right_row) ? $right_row->{acc} : "");
    
    while (1) {
        
        if ( (! $left_row) || (! $right_row)) {
            last;
        }
        
        my $left_acc = $left_row->{acc};
        my $left_median = $left_row->{median_cov};
        my $left_mean = $left_row->{mean_cov}; 
        my $left_stdev = $left_row->{stdev};
            
        my $right_acc = $right_row->{acc};
        my $right_median = $right_row->{median_cov};
        my $right_mean = $right_row->{mean_cov};
        my $right_stdev = $right_row->{stdev};

                
        my $core_acc = $left_acc;
        if ($left_acc =~ /^(\S+)\/\d$/) {
            $core_acc = $1;
        }

        my $right_core = $right_acc;
        if ($right_acc =~ /^(\S+)\/\d$/) {
            $right_core = $1;
        }

        unless ($right_core eq $core_acc) {
            
            if ($sorted_flag) {
                if ($left_acc lt $right_acc) {
                    $left_row = $left_delim_parser->get_row(); # advance left
                    $num_left++ if $left_row;
                }
                else {
                    # advance right
                    $right_row = $right_delim_parser->get_row();
                    $num_right++ if $right_row;
                }
                next;
            }
            else {
                die "Error, core accs are not equivalent: [$core_acc] vs. [$right_core] reads, and --sorted flag wasn't used here.";
            }
        }
        
        
        print join("\t", $core_acc,
                   $left_acc, $left_median, $left_mean, $left_stdev,
                   $right_acc, $right_median, $right_mean, $right_stdev,
                   sprintf("%.1f", ($left_median + $right_median)/2), 
                   sprintf("%.1f", ($left_mean + $right_mean)/2), 
                   sprintf("%.1f", ($left_stdev + $right_stdev)/2)
            ) . "\n";
        
        
        $num_paired++;
        
        ## advance next entries.
        
        $left_row = $left_delim_parser->get_row();
        $right_row = $right_delim_parser->get_row();
        $num_left++ if $left_row;
        $num_right++ if $right_row;
    }

    # count any remaining unpaired entries
    while ($left_delim_parser->get_row()) { $num_left++; }
    while ($right_delim_parser->get_row()) { $num_right++; }
    
    print STDERR "-paired $num_paired reads (left entries: $num_left, right entries: $num_right)\n";
    
    my $read_naming_hint = "Left and right read names must match except for a trailing /1 or /2 (or an Illumina '1:N:...' / '2:N:...' comment), "
        . "and both files must contain the same set of read pairs. Example names seen: left [$example_left_acc], right [$example_right_acc]. "
        . "Reads from SRA should be dumped with the paired reads split into separate files with consistent names "
        . "(eg. fasterq-dump --split-files), and both files should be trimmed/filtered together as pairs. "
        . "You can check your reads with \$TRINITY_HOME/util/misc/validate_fastqs.py "
        . "or run with --no_normalize_reads to skip normalization.\n";
    
    if ($num_paired == 0 && ($num_left > 0 || $num_right > 0)) {
        die "Error, no left/right read pairs could be matched by name during normalization. $read_naming_hint";
    }
    elsif ($num_paired < 0.5 * min($num_left, $num_right)) {
        print STDERR "WARNING: only $num_paired of $num_left left and $num_right right reads could be paired by name. $read_naming_hint";
    }
    
    exit(0);
}
        
