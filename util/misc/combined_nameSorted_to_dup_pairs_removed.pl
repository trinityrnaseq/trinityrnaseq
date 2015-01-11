#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 combined_nameSorted.sam\n\n";

my $combined_nameSorted_sam = $ARGV[0] or die $usage;


main: {

    my $conglom_file = "$combined_nameSorted_sam.__conglom_tmp";
    
    unless (-s $conglom_file) {
        
        open (my $ofh, ">$conglom_file") or die "Error, cannot write to $conglom_file";
        
        my @paired;
        
        my $prev_acc = "";
        
        open (my $fh, $combined_nameSorted_sam) or die $!;
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            if ($x[0] ne $prev_acc) {
                if (scalar(@paired) == 2) {
                    &process_pair($ofh, \@paired);
                }
                @paired = ();
            }
            $prev_acc = $x[0];
            push (@paired, [@x]);
        }
        
        ## get last ones
        if (scalar(@paired) == 2) {
            &process_pair($ofh, \@paired);
        }
        
        close $ofh;
    }

    
    
    my $cmd = "sort -k1,1 -k2,2n -k3,3 -k4,4n $conglom_file > $conglom_file.sorted";
    &process_cmd($cmd) unless (-s "$conglom_file.sorted");

    my $dups_removed_file = "$conglom_file.dups_removed.sam";
    unless (-s $dups_removed_file) {
        open (my $ofh, ">$dups_removed_file") or die $!;
        open (my $fh, "$conglom_file.sorted") or die $!;
        print STDERR "-writing $dups_removed_file\n";
        my $prev_contig_info = "";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $contig_info = join("\t", @x[0..3]);
            
            my $read_A_info = $x[4];
            my $read_B_info = $x[5];
            if ($contig_info ne $prev_contig_info) {
                $read_A_info =~ s/$;/\t/g;
                $read_B_info =~ s/$;/\t/g;

                print $ofh "$read_A_info\n$read_B_info\n";
            }
            $prev_contig_info = $contig_info;
        }
        close $fh;
        close $ofh;
    }
    
    
    
    exit(0);
    

}


####
sub process_pair {
    my ($ofh, $pairs_aref) = @_;
    
    my ($read_A, $read_B) = @$pairs_aref;
    
    if ($read_A->[2] gt $read_B->[2]) {
        ## swap 'em
        ($read_A, $read_B) = ($read_B, $read_A);
    }
    
    my $contig_A = $read_A->[2];
    my $pos_A = $read_A->[3];

    my $contig_B = $read_B->[2];
    my $pos_B = $read_B->[3];
    
    my $read_A_text = join("$;", @$read_A);
    my $read_B_text = join("$;", @$read_B);

    print $ofh join("\t", $contig_A, $pos_A, $contig_B, $pos_B, $read_A_text, $read_B_text) . "\n";
    
    return;
    
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) { 
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}
