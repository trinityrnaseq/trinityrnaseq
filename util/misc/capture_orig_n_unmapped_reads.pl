#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use Fastq_reader;
use SAM_reader;
use SAM_entry;

my $usage = <<_EOUSAGE_;


############################################################################################################
#
#  --aligned_sam <string>      aligned reads in sam format (sorting doesn't matter)
#
#  --sampled_fq_list <string>    comma-delimited list of sampled files (no spaces) 
#                                          eg. "sample_left.fq,sample_right.fq"  or "ssample_single.fq"
#
#  --all_fq_list <string>        same as above, corresponding to the original (complete) set of fastq reads
#
#############################################################################################################


_EOUSAGE_


    ;


my $aligned_sam_file;
my $sampled_fq_list;
my $all_fq_list;


&GetOptions( 'aligned_sam=s' => \$aligned_sam_file,
             'sampled_fq_list=s' => \$sampled_fq_list,
             'all_fq_list=s' => \$all_fq_list,
             );


if (@ARGV) {
    die $usage;
}

unless ($aligned_sam_file && $sampled_fq_list && $all_fq_list) {
    die $usage;
}



main: {

    ## get the list of reads that were in the original sample
    my %sampled_reads;
    foreach my $fq_file (split(/,/, $sampled_fq_list)) {
        $fq_file =~ s/\s+//g;
        
        &get_sampled_read_names(\%sampled_reads, $fq_file);
    }

    my $num_sampled_reads = scalar(keys %sampled_reads);
    print STDERR "$num_sampled_reads  number of sampled reads\n";


    ## identify those reads that are aligned
    my %aligned_reads;
    &parse_aligned_reads(\%aligned_reads, $aligned_sam_file);
    
    my $num_aligned_reads = scalar(keys %aligned_reads);
    print STDERR "$num_aligned_reads number of aligned reads\n";

    
    ## report the unmapped reads and the original ones.
    foreach my $fq_file (split(/,/, $all_fq_list)) {
        $fq_file =~ s/\s//g;
        
        &report_orig_n_unmapped(\%sampled_reads, \%aligned_reads, $fq_file);
    }


    exit(0);
    
}


####
sub get_sampled_read_names {
    my ($sampled_reads_href, $fq_file) = @_;

    my $x = 0;
    my $fq_parser = new Fastq_reader($fq_file);
    while (my $record = $fq_parser->next()) {
        my $full_read_name = $record->get_full_read_name();
        #print "SAMPLED: [$full_read_name]\n";
        $sampled_reads_href->{$full_read_name}++;
    
        $x++;
        #if ($x > 5) { last; }
    }

    return;
}


####
sub parse_aligned_reads {
    my ($aligned_reads_href, $aligned_sam_file) = @_;
    
    my $sam_reader = new SAM_reader($aligned_sam_file);
    
    my $num_reads_aligned = 0;

    while (my $sam_entry = $sam_reader->get_next()) {

        if ($sam_entry->is_query_unmapped()) { next; }

        my $read_name = $sam_entry->reconstruct_full_read_name();
        
        #print "SAM: $read_name\n";

        $aligned_reads_href->{$read_name}++;
        
        $num_reads_aligned++;

        
        #if ($num_reads_aligned > 5) { last; }
    }

    #print "$aligned_sam_file contains $num_reads_aligned aligned reads\n\n";

    return;
}

####
sub report_orig_n_unmapped {
    my ($sampled_reads_href, $aligned_reads_href, $fq_file) = @_;
    
    my $fq_parser = new Fastq_reader($fq_file);

    my $num_reads_output = 0;
    my $num_reads_skipped = 0;
    
    

    while (my $record = $fq_parser->next()) {
        
        my $read_name = $record->get_full_read_name();
        
        #print "FQ: [$read_name]\n";
        
        my $print_record_flag = $sampled_reads_href->{$read_name} || 0;
        

        #print "sampled: $read_name = $print_record_flag\n";
        
        #$print_record_flag = 0;
        
        unless ($print_record_flag) {

            if ($read_name =~ /^(\S+)\/([12])$/) {

                ## check to see if both read pairs align. If not, report both.

                my $core = $1;
                my $pair_val = $2;

                my $other_val = ($pair_val == 1) ? 2 : 1;

                my $other_read_name = join("/", $core, $other_val);

                my $got_aligned_read = $aligned_reads_href->{$read_name} || 0;
                my $got_other_aligned_read = $aligned_reads_href->{$other_read_name} || 0;

                #print "aligned_read: $read_name = $got_aligned_read\n";
                #print "other_read: $other_read_name = $got_other_aligned_read\n";
                

                unless ($aligned_reads_href->{$read_name} && $aligned_reads_href->{$other_read_name}) {
                    $print_record_flag = 1;
                }
            }
            else {
                # single read
                
                my $got_single = $aligned_reads_href->{$read_name} || 0;
                #print "single: $read_name = $got_single\n";
                
                unless ($aligned_reads_href->{$read_name}) {
                    $print_record_flag = 1;
                }
            }
        }

        if ($print_record_flag) {
            my $fastq_text = $record->get_fastq_record();
            print $fastq_text;
            
            $num_reads_output++;
        }
        else {
            $num_reads_skipped++;
        }

        #if ($num_reads_output > 5) { last; }
        
    }
    
    print STDERR "\n$fq_file: Total number of reads output: $num_reads_output, num skipped: $num_reads_skipped = " 
        . sprintf("%.2f", $num_reads_skipped/($num_reads_output+$num_reads_output)*100) . "\n";
    
    return;
}



        
