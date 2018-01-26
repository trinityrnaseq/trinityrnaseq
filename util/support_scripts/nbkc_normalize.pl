#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $pair_stats_file;
my $max_cov;
my $min_cov;
my $max_pct_stdev;


my $usage = <<__EOUSAGE__;

#############################################################
#
# Required:
#
#  --stats_file <string>     : pairs.stats.sorted
#
#  --max_cov <int>           : maximum coverage
#
#  --min_cov <int>           : minimum coverage
#
#  --max_pct_stdev <int>     : maximum pct stdev
#
#
#############################################################

__EOUSAGE__

    ;



my $help_flag;

&GetOptions ( 'h' => \$help_flag,

              'stats_file=s' => \$pair_stats_file,
              'max_cov=i' => \$max_cov,
              'min_cov=i' => \$min_cov,
              'max_pct_stdev=i' => \$max_pct_stdev,

              
);


unless ($pair_stats_file &&
        $max_cov &&
        defined($min_cov) &&
        $max_pct_stdev) {
    
    die $usage;

}



main: {

    my $count_no_cov = 0;
    my $count_aberrant_and_discarded = 0;
    my $count_selected = 0;
    my $count_total = 0;

    my $count_below_min_cov = 0;
    
    open (my $fh, $pair_stats_file) or die $!;
    while (<$fh>) {
        
        $count_total++;

        chomp;
        my $line = $_;

        my ($med_cov, $avg_cov, $stdev, $pct_dev, $core_acc) = split(/\t/);
        
        $core_acc =~ s|/[12]$||;
        
        if ($med_cov < 1) { 
            
            $count_no_cov++;   ## this shouldn't happen in modern versions of this process. (bh 10-2013)
            next; 
        }

        if ($med_cov < $min_cov) {
            $count_below_min_cov++;
            next;
        }
        
        if ($pct_dev > $max_pct_stdev) { 
            
            # print STDERR "// excluding $_ as aberrant\n";
            
            $count_aberrant_and_discarded++;
            next; 
        }
                                
        if (rand(1) <= $max_cov/$med_cov) {
            print "$core_acc\n";
            $count_selected++;
        }
    }
    close $fh;
    
    print STDERR "$count_selected / $count_total = " . sprintf("%.2f", $count_selected/$count_total*100) . "% reads selected during normalization.\n";
    print STDERR "$count_aberrant_and_discarded / $count_total = " . sprintf("%.2f", $count_aberrant_and_discarded/$count_total*100) . "% reads discarded as likely aberrant based on coverage profiles.\n";
    print STDERR "$count_no_cov / $count_total = " . sprintf("%.2f", $count_no_cov/$count_total*100) . "% reads missing kmer coverage (N chars included?).\n";
    print STDERR "$count_below_min_cov / $count_total = " . sprintf("%.2f", $count_below_min_cov/$count_total*100) . "% reads discarded as below minimum coverage threshold=$min_cov\n";
    
    exit(0);
}


