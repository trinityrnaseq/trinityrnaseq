#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");

use DelimParser;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $pair_stats_file;
my $max_cov;
my $min_cov=1;
my $max_CV;


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
#  --max_CV <int>            : maximum coeff. var.
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
              'max_CV=i' => \$max_CV,

              
);


unless ($pair_stats_file &&
        $max_cov &&
        defined($min_cov) &&
        defined($max_CV)) {
    
    die $usage;

}



main: {

    my $count_aberrant_and_discarded = 0;
    my $count_selected = 0;
    my $count_total = 0;

    my $count_below_min_cov = 0;

    open (my $fh, $pair_stats_file) or die $!;

    my $delim_parser = new DelimParser::Reader($fh, "\t");
    
    while (my $row = $delim_parser->get_row()) {
        
        $count_total++;

        my $core_acc = $row->{acc};
        
        $core_acc =~ s|/[12]$||;

        my $med_cov = ($row->{median_cov});

        if ($med_cov < $min_cov) {
            $count_below_min_cov++;
            next;
        }

        my $sd = $row->{stdev};
        my $u = $row->{mean_cov};

        if ($u <= 0) {
            $count_aberrant_and_discarded++;
            next; 
        }   
        
        my $cv = $sd/$u;
        
        if ($cv > $max_CV) {
            $count_aberrant_and_discarded++;
            next; 
        }
                                
        if (rand(1) <= $max_cov/$med_cov) {
            print "$core_acc\n";
            $count_selected++;
        }
    }
    close $fh;

    unless ($count_total) {
        die "Error, no reads made it to the normalization process...  ";
    }
    
    print STDERR "$count_selected / $count_total = " . sprintf("%.2f", $count_selected/$count_total*100) . "% reads selected during normalization.\n";
    print STDERR "$count_aberrant_and_discarded / $count_total = " . sprintf("%.2f", $count_aberrant_and_discarded/$count_total*100) . "% reads discarded as likely aberrant based on coverage profiles.\n";
    print STDERR "$count_below_min_cov / $count_total = " . sprintf("%.2f", $count_below_min_cov/$count_total*100) . "% reads discarded as below minimum coverage threshold=$min_cov\n";
    
    exit(0);
}


