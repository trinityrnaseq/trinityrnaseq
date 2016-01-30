#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fastq_reader;

my $DEBUG = 0;


my $usage = "usage: $0 left.fq right.fq\n\n";

my $left_fq = $ARGV[0] or die $usage;
my $right_fq = $ARGV[1] or die $usage;

open (my $left_PP_ofh, ">$left_fq.P.fq") or die $!;
open (my $left_UP_ofh, ">$left_fq.U.fq") or die $!;

open (my $right_PP_ofh, ">$right_fq.P.fq") or die $!;
open (my $right_UP_ofh, ">$right_fq.U.fq") or die $!;

my $ok_counter = 0;
my $left_orphan_counter = 0;
my $right_orphan_counter = 0;

main: {

    my $left_fq_reader = new Fastq_reader($left_fq);
    my $right_fq_reader = new Fastq_reader($right_fq);

        
        
    my ($left_fq_record, $right_fq_record);
    
    my @left_entries;
    my @right_entries;
    
    my %core_counter;

    do {
        

        my $num_left_stored = scalar(@left_entries);
        my $num_right_stored = scalar(@right_entries);


        if ($DEBUG) {
        
            my %seen;
            
            foreach my $left_entry (@left_entries) {
                print STDERR "L " . $left_entry->get_full_read_name() . "\n" if $DEBUG;
                my $core_acc = $left_entry->get_core_read_name();
                $seen{$core_acc}++;
            }
            print STDERR "\n" if $DEBUG;
            my $found_hit = 0;
            foreach my $right_entry (@right_entries) {
                print STDERR "R " . $right_entry->get_full_read_name() . "\n" if $DEBUG;
                my $core_acc = $right_entry->get_core_read_name();
                my $count = ++$seen{$core_acc};
                if ($count == 2) {
                    print STDERR " ***** \n" if $DEBUG;
                    $found_hit++;
                }
            }
            print STDERR "\n\n" if $DEBUG;
            
            if ($found_hit) {
                die " reads must be jumbled";
            }
        }
        
        my $MAX_ORPHAN_STORE = 100;
        #if ($num_left_stored > $MAX_ORPHAN_STORE &&  $num_right_stored > $MAX_ORPHAN_STORE) { die; }
        
        
        if ($ok_counter % 1000 == 0) {
            print STDERR "\r[$ok_counter pairs_written, left_orphans_written: $left_orphan_counter, right_orphans_written: $right_orphan_counter]    ";
            print STDERR "[Left cache:$num_left_stored, Right cache:$num_right_stored]    ";
        }
        
        $left_fq_record = $left_fq_reader->next();
        push (@left_entries, $left_fq_record);
        
        my $left_core_acc = $left_fq_record->get_core_read_name();
        my $count = ++$core_counter{$left_core_acc};
        
        if ($count == 2) {
            &dump_pairs($left_core_acc, \@left_entries, \@right_entries, \%core_counter);
        }
        
        $right_fq_record = $right_fq_reader->next();
        push (@right_entries, $right_fq_record);

        my $right_core_acc = $right_fq_record->get_core_read_name();
        $count = ++$core_counter{$right_core_acc};
        
        if ($count == 2) {
            &dump_pairs($right_core_acc, \@left_entries, \@right_entries, \%core_counter);
        }


        #print STDERR Dumper(\%core_counter);
        
        


    } while ($left_fq_record && $right_fq_record);
    
    while (@left_entries) {
        
        $left_fq_record = shift @left_entries;
        
        print $left_UP_ofh $left_fq_record->get_fastq_record();
        
        $left_orphan_counter++;
    }
    
    while (@right_entries) {
        
        $right_fq_record = shift @right_entries;
        
        print $right_UP_ofh $right_fq_record->get_fastq_record();
    
        $right_orphan_counter++;
        
    }

    print STDERR "\r[$ok_counter pairs_written, left_orphans_written: $left_orphan_counter, right_orphans_written: $right_orphan_counter]\n\nDone.\n\n";    
    
    exit(0);
}
        
####
sub dump_pairs {
    my ($acc, $left_entries_aref, $right_entries_aref, $core_counter_href) = @_;

    if ($left_entries_aref->[ $#$left_entries_aref ]->get_core_read_name() eq $acc) {
        
        my $record = pop @$left_entries_aref;
        print $left_PP_ofh $record->get_fastq_record();
        
        ## write earlier stored records as unpaired entries
        while ($record = shift @$left_entries_aref) {
            
            my $core_acc = $record->get_core_read_name();
            delete $core_counter_href->{$core_acc};
            
            print $left_UP_ofh $record->get_fastq_record();
            $left_orphan_counter++;
        }
        
        # process right records
        while ($record = shift @$right_entries_aref) {
            
            my $core_acc = $record->get_core_read_name();
            if ($core_acc eq $acc) {
                print $right_PP_ofh $record->get_fastq_record();
                last; # retain any remaining entries on the stack 
                
            }
            else {
                
                my $core_acc = $record->get_core_read_name();
                delete $core_counter_href->{$core_acc};
                print $right_UP_ofh $record->get_fastq_record();
            
                $right_orphan_counter++;
            }
            
        }
    }
    elsif ($right_entries_aref->[ $#$right_entries_aref ]->get_core_read_name() eq $acc) {
        
        my $record = pop @$right_entries_aref;
        print $right_PP_ofh $record->get_fastq_record();
        
        ## write earlier stored records as unpaired entries
        while ($record = shift @$right_entries_aref) {
            
            my $core_acc = $record->get_core_read_name();
            delete $core_counter_href->{$core_acc};
            
            print $right_UP_ofh $record->get_fastq_record();
            
            $right_orphan_counter++;
        }
        
        # process left records
        while ($record = shift @$left_entries_aref) {
            
            my $core_acc = $record->get_core_read_name();
            if ($core_acc eq $acc) {
                print $left_PP_ofh $record->get_fastq_record();
                last; # retain any remaining entries on the stack 
                
            }
            else {
                
                my $core_acc = $record->get_core_read_name();
                delete $core_counter_href->{$core_acc};
                print $left_UP_ofh $record->get_fastq_record();
                
                $left_orphan_counter++;
            }
            
        }
    }
    
    delete $core_counter_href->{$acc};

    print STDERR "\n\nOK: $acc\n" if $DEBUG;
    
    $ok_counter += 2;
    
    return;
}

