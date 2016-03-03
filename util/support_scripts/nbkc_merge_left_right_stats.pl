#!/usr/bin/env perl

use strict;
use warnings;

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
    
    my $left_text = <$left_fh>;
    my $right_text = <$right_fh>;

    while (1) {

        #print STDERR "$left_text\n$right_text\n\n";

                
        if ( (! $left_text) || (! $right_text)) {
            last;
        }


                
        chomp $left_text;
        chomp $right_text;

        my ($left_median, $left_avg, $left_stdev, $left_pct, $left_acc, $left_cov) = split(/\t/, $left_text);
        my ($right_median, $right_avg, $right_stdev, $right_pct, $right_acc, $right_cov) = split(/\t/, $right_text);
        
        my $core_acc = $left_acc;
        if ($left_acc =~ /^(\S+)\/\d$/) {
            $core_acc = $1;
        }

        my $right_core = $right_acc;
        if ($right_acc =~ /^(\S+)\/\d$/) {
            $right_core = $1;
        }

        #print STDERR "$core_acc\t$right_core";
        #if ($core_acc eq $right_core) { print STDERR "\tYES\n";} else { print STDERR "\n"; }
        
        unless ($right_core eq $core_acc) {
            
            if ($sorted_flag) {
                if ($left_acc lt $right_acc) {
                    $left_text = <$left_fh>; # advance left
                }
                else {
                    # advance right
                    $right_text = <$right_fh>;
                }
                next;
            }
            else {
                die "Error, core accs are not equivalent: [$core_acc] vs. [$right_core] reads, and --sorted flag wasn't used here.";
            }
        }
        
        my $median_pair_cov = int( ($left_median+$right_median)/2 + 0.5);
        my $avg_pair_cov = int ( ($left_avg + $right_avg)/2 + 0.5);
        my $avg_stdev = int( ($left_stdev + $right_stdev)/2 + 0.5);
        
        if ($median_pair_cov != 0) { 
            ## ignore bad reads.
            
            my $avg_pair_pct = int( ($left_pct + $right_pct) / 2 + 0.5);
            
            print join("\t", $median_pair_cov, $avg_pair_cov, $avg_stdev, $avg_pair_pct, $core_acc) . "\n";
        }
    

        ## advance next entries.
        
        $left_text = <$left_fh>;
        $right_text = <$right_fh>;
        
    }
    
    exit(0);
}
        
