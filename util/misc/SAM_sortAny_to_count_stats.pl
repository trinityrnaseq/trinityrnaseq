#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\n\nusage: aligned_reads.sam [debug]\n\n"
    . "Note: uses the bitflag settings for counting entries. Proper-pairing supercedes other settings.\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $DEBUG = $ARGV[1] || 0;


my $DEBUG_OFH;
if ($DEBUG) {
    open ($DEBUG_OFH, ">_debug.frag_classes") or die $!;
}

=notes

Entirely works off flag settings.

Proper pairs trump individual left/right alignments.

Improper pairs: left and right read alignments exist anywhere but not flagged as proper pairs.

=cut

main: {

    my $prev_read_name = "";
    my $prev_scaff_name = "";

    my %counts;

    my $count = 0;

    my $sam_reader = new SAM_reader($sam_file);
    while ($sam_reader->has_next()) {
        
        $count++;

        my $read = $sam_reader->get_next();
        
        if ($read->is_query_unmapped()) { next; } # ignore unaligned read entries.

        my $scaff_name = $read->get_scaffold_name();
        my $core_read_name = $read->get_core_read_name();
        
        if ($read->is_paired()) {
            
            if ($read->is_proper_pair()) {
                ## erase any single entries.
                if (exists $counts{ $core_read_name . "::L" }) {
                    delete $counts{ $core_read_name . "::L" };
                }
                if (exists $counts{ $core_read_name . "::R" } ) {
                    delete $counts{ $core_read_name . "::R" };
                }
                $counts{ $core_read_name . "::PP" } = 1;
            }
            else {
                # not propper pair
                unless (exists $counts{ $core_read_name . "::PP" }) {
                    if ($read->is_first_in_pair()) {
                        $counts{ $core_read_name . "::L" } = 1;
                    }
                    elsif ($read->is_second_in_pair()) {
                        $counts{ $core_read_name . "::R" } = 1;
                    }
                }
            }
        }
        else {
            $counts{ $core_read_name . "::S" } = 1;
        }
        
        if ($count % 100000 == 0) {
            my $count_print = $count;
            $count_print=~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
            print STDERR "\r[$count_print]   ";
        }
    }
    print STDERR "\n\n";
    
    ## identify improper pairs
    my @reads = keys %counts;
    foreach my $read (@reads) {
        if ($read =~ /::L$/) {
            my $core_name = $read;
            $core_name =~ s/::L$//;
            
            if (exists($counts{ $core_name . "::R" }) ) {
                ## count as improper pair instead
                $counts{ $core_name . "::IP" } = 1;
                delete $counts{$read};
                delete $counts{ $core_name . "::R" };
                
            }
        }
    }


    ## sum counts, generate summary
    
    my $count_PP = 0;
    my $count_L = 0;
    my $count_R = 0;
    my $count_S = 0;
    my $count_IP = 0; # improper pairs

    my $total = 0;

    foreach my $read (keys %counts) {
        
        my @x = split(/::/, $read);
        my $class = pop @x;
        my $core_read_name = join("::", @x);
        
        if ($class eq "PP") {
            $count_PP += 2;
            $total += 2;
        }
        elsif ($class eq "IP") {
            $count_IP += 2;
            $total += 2;

        }
        elsif ($class eq "L") {
            $count_L++;
            $total++;

        }
        elsif ($class eq "R") {
            $count_R++;
            $total++;

        }
        elsif ($class eq "S") {
            $count_S++;
            $total++;

        }
    
        if ($DEBUG) {
            print $DEBUG_OFH join("\t", $core_read_name, $class) . "\n";
                            
        }
        
    }
    
    print "Proper_pair:\t$count_PP\t" . sprintf("%.3f%%", $count_PP/$total*100) . "\n"
        . "Improper_pair:\t$count_IP\t" . sprintf("%.3f%%", $count_IP/$total*100) . "\n"
        . "Left-only:\t$count_L\t" . sprintf("%.3f%%", $count_L/$total*100) . "\n"
        . "Right-only:\t$count_R\t" . sprintf("%.3f%%", $count_R/$total*100) . "\n"
        . "Single_end:\t$count_S\t" . sprintf("%.3f%%", $count_S/$total*100) . "\n\n";

    print "Total aligned frags: $total\n\n";



    close $DEBUG_OFH if $DEBUG;

    exit(0);
    

}


    

        
    
