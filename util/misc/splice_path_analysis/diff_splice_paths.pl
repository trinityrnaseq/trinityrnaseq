#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 splicePaths_A  splicePaths_B\n\n";

my $fileA = $ARGV[0] or die $usage;
my $fileB = $ARGV[1] or die $usage;


main: {

    my %combosA = &get_combos($fileA);
    my %combosB = &get_combos($fileB);

    my %introns_to_combos;
    &index_combos_by_intron(\%combosA, \%introns_to_combos);
    &index_combos_by_intron(\%combosB, \%introns_to_combos);

    my %seen_intron;

    foreach my $intron (keys %introns_to_combos) {
        
        if ($seen_intron{$intron}) { next; }
                
        my %paths;
        &get_all_connected_paths($intron, \%introns_to_combos, \%seen_intron, \%paths);
        
        

        my @paths_A_not_B;
        my @paths_B_not_A;
        my @paths_both_A_and_B;
        foreach my $path (keys %paths) {
            if ($combosA{$path} && $combosB{$path}) {
                push (@paths_both_A_and_B, $path);
            }
            elsif ($combosA{$path}) {
                push (@paths_A_not_B, $path);
            }
            elsif ($combosB{$path}) {
                push (@paths_B_not_A, $path);
            }
        }


        ## filter out those that are subpaths of others.
        @paths_A_not_B = &remove_subpaths(\@paths_A_not_B, [@paths_B_not_A]);
        @paths_B_not_A = &remove_subpaths(\@paths_B_not_A, [@paths_A_not_B]);
        

        print join("\t", $intron, 
                   scalar(@paths_both_A_and_B),
                   scalar(@paths_A_not_B),
                   scalar(@paths_B_not_A)) . "\n";
    }
    
    
    
    exit(0);
    
}

####
sub remove_subpaths {
    my ($paths_query_aref, $paths_to_examine_as_containing_subpaths_aref) = @_;

    my @ok;

    foreach my $path (@$paths_query_aref) {
       
        my ($chr, $rest_path) = split(/:/, $path);

        my $found_as_subpath = 0;

        foreach my $other_path (@$paths_to_examine_as_containing_subpaths_aref) {

            if ($other_path =~ /$rest_path/) {
                $found_as_subpath = 1;
                last;
            }
        }
        if (! $found_as_subpath) {
            push (@ok, $path);
        }
    }

    return(@ok);
}


####
sub get_combos {
    my ($file) = @_;

    my %combos;

    open (my $fh, $file) or die $!;
    while (<$fh>) {
        chomp;
        
        my @x = split(/\t/);
        
        if ($x[0] =~ /COMBO/) { 
            $combos{$x[1]} = 1;
        }

    }

    close $fh;

    return(%combos);
}


####
sub index_combos_by_intron {
    my ($combos_href, $introns_to_combos_href) = @_;
    
    foreach my $combo (keys %$combos_href) {
        my @introns = &get_introns_from_path($combo);        
        foreach my $intron (@introns) {
            push (@{$introns_to_combos_href->{$intron}}, $combo);
        }
    }
    
    return;
}


sub get_introns_from_path {
    my ($path) = @_;
    my ($chr, @introns) = split(/_+/, $path);
        
    foreach my $intron (@introns) {
        $intron = "$chr$intron";
    }
    return(@introns);
}


sub get_all_connected_paths {
    my ($intron, $introns_to_combos_href, $seen_href, $paths_href) = @_;
    
    my @paths = @{$introns_to_combos_href->{$intron}};
    $seen_href->{$intron} = 1;
    foreach my $path (@paths) {
        if (! exists $paths_href->{$path}) {
            $paths_href->{$path} = 1;
            my @introns = &get_introns_from_path($path);
            foreach my $intron (@introns) {
                if ($seen_href->{$intron}) { next; }
                &get_all_connected_paths($intron, $introns_to_combos_href, $seen_href, $paths_href);
            }
        }
    }

    return;
}
