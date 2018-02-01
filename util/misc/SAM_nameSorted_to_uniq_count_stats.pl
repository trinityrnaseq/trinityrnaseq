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

    my %counts = ( 'UPP' => 0, # unique proper pairs
                   'MPP' => 0, # multi-mapped proper pairs
                   
                   'IP' => 0, # improper pairs
                   
                   'UL' => 0, # left unique reads
                   'ML' => 0, # left multi-mapped reads
                   
                   'UR' => 0, # right unique reads
                   'MR' => 0, # right multi-mapped reads

                   'US' => 0, # single unique-reads
                   'MS' => 0, # single multi-mapped reads
        );
         

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

=try_mirror_bowtie2_summary

30575 reads; of these:
  30575 (100.00%) were paired; of these:
    2577 (8.43%) aligned concordantly 0 times
    5766 (18.86%) aligned concordantly exactly 1 time
    22232 (72.71%) aligned concordantly >1 times
    ----
    2577 pairs aligned concordantly 0 times; of these:
      133 (5.16%) aligned discordantly 1 time
    ----
    2444 pairs aligned 0 times concordantly or discordantly; of these:
      4888 mates make up the pairs; of these:
        3588 (73.40%) aligned 0 times
        322 (6.59%) aligned exactly 1 time
        978 (20.01%) aligned >1 times
94.13% overall alignment rate

=cut

    print "Stats for aligned rna-seq fragments (note, not counting those frags where neither left/right read aligned)\n";
    print "\n\n$sum_reads aligned fragments; of these:\n";
    print "  " . ($sum_reads - $counts{US} - $counts{MS}) . " were paired; of these:\n";
    print "    " . ($counts{UL} + $counts{ML} + $counts{UR} + $counts{MR} + $counts{IP}) . " aligned concordantly 0 times\n";
    print "    " . $counts{UPP} . " aligned concordantly exactly 1 time\n";
    print "    " . $counts{MPP} . " aligned concordantly >1 times\n";
    print "    ----\n";
    print "    " . ($counts{UL} + $counts{ML} + $counts{UR} + $counts{MR} + $counts{IP}) . " pairs aligned concordantly 0 times; of these:\n";
    print "    " . $counts{IP} . " aligned as improper pairs\n";
    print "    " . ($counts{UL} + $counts{ML} + $counts{UR} + $counts{MR}) . " pairs had only one fragment end align to one or more contigs; of these:\n";
    print "       " . ($counts{UL} + $counts{ML}) . " fragments had only the left /1 read aligned; of these:\n";
    print "            " . $counts{UL} . " left reads mapped uniquely\n";
    print "            " . $counts{ML} . " left reads mapped >1 times\n";
    print "       " . ($counts{UR} + $counts{MR}) . " fragments had only the right /2 read aligned; of these:\n";
    print "            " . $counts{UR} . " right reads mapped uniquely\n";
    print "            " . $counts{MR} . " right reads mapped >1 times\n";
    print "Overall,  " . sprintf("%.2f", ($counts{UPP} + $counts{MPP})/$sum_reads * 100) . "% of aligned fragments aligned as proper pairs\n\n\n";
    
    if ($counts{US} + $counts{MS}) {
        print "  " . ($counts{US} + $counts{MS}) . " were single-end reads; of these:\n";
        print "      " . $counts{US} . " single-end reads mapped uniquely\n";
        print "      " . $counts{MS} . " single-end reads mapped >1 times\n";
    }
        
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
    my @proper_pairs;
    
    
  pair_search:
    foreach my $left_read (@left_reads) {

        unless ($left_read->is_proper_pair()) { next; }
        
        
        my $aligned_pos = $left_read->get_aligned_position();
                
        foreach my $right_read (@right_reads) {

            unless ($right_read->is_proper_pair()) { next; }
            
            if ($left_read->get_scaffold_name() eq $right_read->get_scaffold_name()
                &&
                $right_read->get_mate_scaffold_position() == $aligned_pos) {
                
                push (@proper_pairs, [$left_read, $right_read]);
                
                last;
            }
        }
    }
    

    my $class = "";
    
    if (@proper_pairs) {

        if (scalar(@proper_pairs) == 1) {
            $class = "UPP";
        }
        else {
            # multi ampped proper pairs
            $class = "MPP";
        }
    }
    elsif ($got_left_read && $got_right_read) {
        $class = "IP";
    }
    elsif ($got_left_read) {

        if (scalar(@left_reads) == 1) {
            $class = "UL";
        }
        else {
            $class = "ML";
        }
    }
    elsif ($got_right_read) {
        
        if (scalar(@right_reads) == 1) {
            $class = "UR";
        }
        else {
            $class = "MR";
        }
    }
    else {
        if (scalar(@reads) == 1) {
            $class = "US";
        }
        else {
            $class = "MS";
        }
    }

    $counts_href->{$class}++;
    
    if ($DEBUG) {
        my $core_read_name = $reads[0]->get_core_read_name();
        
        print $DEBUG_OFH join("\t", $core_read_name, $class) . "\n";
    }

    return;
}


    

        
    
