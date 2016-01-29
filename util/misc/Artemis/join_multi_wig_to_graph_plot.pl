#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use WigParser;


my $usage = "usage: $0 wigA [ wigB ... ]\n\n";
my @wig_files = @ARGV;


unless (@wig_files) {
    die $usage;
}


main: {

    my @data_arefs;
    
    my $max_pos = 0;
    foreach my $wig_file (@wig_files) {
        
        my %scaff_to_data = &WigParser::parse_wig($wig_file);
        
        my @scaffs = keys %scaff_to_data;
        if (scalar(@scaffs) != 1) {
            die "Error, only a single scaffold can be leveraged by Artemis, and " . scalar(@scaffs) . " scaffs worth of data are found in file: $wig_file:\n@scaffs";
        }

        my $scaff = $scaffs[0];
        
        my $data_aref = $scaff_to_data{$scaff};
        
        push (@data_arefs, $data_aref);
        
        my $max = $#$data_aref;
        if ($max > $max_pos) {
            $max_pos = $max;
        }


    }

    ## output results

    

    
    for (my $i = 1; $i <= $max_pos; $i++) {
        
        my $printed_flag = 0;
        foreach my $data_aref (@data_arefs) {
            
            my $val = $data_aref->[$i] || 0;
            if ($printed_flag) {
                print "\t";
            }

            print "$val";
            $printed_flag = 1;
            
        }
        print "\n";
    }
                
    
    
    exit(0);
    


    
}
