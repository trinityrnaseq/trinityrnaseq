#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "\n\nusage: name_sorted_paired_reads.sam [debug]\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $DEBUG = $ARGV[1] || 0;

if ($sam_file =~ /coord/i) {
    die "Your filename contains the term 'coord' which makes me nervous.  This needs to be a name-sorted SAM file. Please rename the file or use the correct version here\n.";
}


my $DEBUG_OFH;
if ($DEBUG) {
    open ($DEBUG_OFH, ">_debug.nameSorted_frag_classes") or die $!;
}


=notes

Gathers all reads having the same core read name (both left and right reads of PE frags).

Separates the left and right reads into corresponding lists.

Examines pairwise comparisons between left and right reads, check for same scaffold and that the right reads mate aligned coordinate matches up with that of the left reads aligned coordinate.

If not properly paired:
    have both left and right reads mapped anywhere:  improper pair
Otherwise, left-only or right-only counts.


=cut





main: {


    my $prev_read_name = "";
    my $prev_scaff_name = "";

    my @reads;

    my %counts;

    my $line_counter = 0;

    my $sam_reader = new SAM_reader($sam_file);
    while ($sam_reader->has_next()) {
        
        $line_counter++;

        my $read = $sam_reader->get_next();
        
        if ($read->is_query_unmapped()) { next; }


        my $scaff_name = $read->get_scaffold_name();
        my $core_read_name = $read->get_core_read_name();
        
        if ($core_read_name ne $prev_read_name) {

            if (@reads) {
                &process_pairs(\@reads, \%counts);
                @reads = ();
            }
        }
        
        push (@reads, $read);
        

        $prev_read_name = $core_read_name;
        $prev_scaff_name = $scaff_name;


        if ($line_counter % 100000 == 0) {
            my $count_print = $line_counter;
            $count_print=~ s/(\d)(?=(\d{3})+(\D|$))/$1\,/g;
            print STDERR "\r[$count_print]  lines read ";
        }
            
        
    }

    print STDERR "\n\n";
        
    &process_pairs(\@reads, \%counts) if @reads;
    

    my $sum_reads = 0;
    foreach my $count (values %counts) {
        $sum_reads += $count;
    }

    print "\n#read_type\tcount\tpct\n";
    foreach my $count_type (reverse sort {$counts{$a}<=>$counts{$b}} keys %counts) {
        my $count = $counts{$count_type};
        print "$count_type\t$count\t" . sprintf("%.2f", $count/$sum_reads*100) . "\n";
    }
    print "\n";
    print "Total aligned reads: $sum_reads\n\n";


    close $DEBUG_OFH if $DEBUG;
    
    exit(0);
}

####
sub process_pairs {
    my ($reads_aref, $counts_href) = @_;
    
    my @reads = @$reads_aref;
    
    my @left_reads;
    my @right_reads;


    foreach my $read (@reads) {
        if ($read->is_first_in_pair()) {
            push (@left_reads, $read);
        }
        elsif ($read->is_second_in_pair()) {
            push (@right_reads, $read);
        }
    }

    my $got_left_read = (@left_reads) ? 1 : 0;
    my $got_right_read = (@right_reads) ? 1: 0;
    

    ## check to see if we have proper pairs:
    my $got_proper_pair = 0;
    
  pair_search:
    foreach my $left_read (@left_reads) {
        
        my $aligned_pos = $left_read->get_aligned_position();
        
        
        foreach my $right_read (@right_reads) {
            
            if ($left_read->get_scaffold_name() eq $right_read->get_scaffold_name()
                &&
                $right_read->get_mate_scaffold_position() == $aligned_pos) {
                
                $got_proper_pair = 1;
                
                last pair_search;
            }
        }
    }


    my $class = "";
    
    if ($got_proper_pair) {
        $counts_href->{proper_pairs} += 2;
        $class = "PP";
        
    }
    elsif ($got_left_read && $got_right_read) {
        $counts_href->{improper_pairs} += 2;
        $class = "IP";
    }
    elsif ($got_left_read) {
        $counts_href->{left_only}++;
        $class = "L";
    }
    elsif ($got_right_read) {
        $counts_href->{right_only}++;
        $class = "R";
    }
    else {
        $counts_href->{single}++;
        $class = "S";
    }
    
    if ($DEBUG) {
        my $core_read_name = $reads[0]->get_core_read_name();
        
        print $DEBUG_OFH join("\t", $core_read_name, $class) . "\n";
    }

    return;
}


    

        
    
