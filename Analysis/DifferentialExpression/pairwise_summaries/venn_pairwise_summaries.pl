#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

#############################################################################################

usage: $0 [--top_ranks <int>] pw1 pw2 ...

##############################################################################################
#
# --top_ranks <int>          restrict to number of top ranked entries per pairwise condition.
#
# --minFC <float>            minimum fold change in expression
#
##############################################################################################


__EOUSAGE__

    ;



my $help_flag;
my $top_ranks = -1;
my $minFC = 0;

&GetOptions ( 'h' => \$help_flag,

              'top_ranks=i' => \$top_ranks,
              'minFC=f' => \$minFC,


              );


my $minLogFC = 0;
if ($minFC) {
    $minLogFC = log($minFC)/log(2);
}


my @summaries = @ARGV or die $usage;

unless (scalar @summaries > 1) {
    die $usage;
}




main: {
    
    
    if ($top_ranks > 0) {
        print STDERR "NOTE: only $top_ranks top ranks from each will be reported.\n\n";
    }
    
    my %gene_samples_to_methods;
    my %gene_samples_to_ranks;

    foreach my $summary (@summaries) {
        
        my %pw_counter;
        open (my $fh, $summary) or die $!;
        while (<$fh>) {
            chomp;
            my ($acc, $sampleA, $sampleB, $exprA, $exprB, $logFC, $pval) = split(/\t/);
            
            if (abs($logFC) < $minLogFC) { next; }
            
            my $feature_sample_key = join("$;", $acc, sort ($sampleA, $sampleB));
            
            my $sample_key = join("$;", sort($sampleA, $sampleB));
            
            $pw_counter{$sample_key}++;
            #if ($top_ranks > 0 && $pw_counter{$sample_key} > $top_ranks) { next; }
            
            #print STDERR "$pw_counter{$sample_key}\n";
            
            push (@{$gene_samples_to_methods{$feature_sample_key}}, $summary);
            
            $gene_samples_to_ranks{$feature_sample_key}->{$summary} = $pw_counter{$sample_key};
            
        }
        close $fh;
    }

    ## output
    
    foreach my $feature_samples (keys %gene_samples_to_methods) {

        my @methods = sort (@{$gene_samples_to_methods{$feature_samples}});
        my @rankings;
        my $report_entry = ($top_ranks > 0) ? 0 : 1;
        foreach my $method (@methods) {
            my $rank = $gene_samples_to_ranks{$feature_samples}->{$method};
            push (@rankings, $rank);

            if ($top_ranks > 0 && $rank <= $top_ranks) {
                $report_entry = 1;
            }
            
        }
        
        if ($report_entry) {
            my ($feature, $sampleA, $sampleB) = split(/$;/, $feature_samples);
            
            print join("\t", $feature, $sampleA, $sampleB, join(",", @methods), join(",", @rankings) ) . "\n";
        }
    }
    
    exit(0);
}


