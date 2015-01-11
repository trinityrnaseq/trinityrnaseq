#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::Descriptive;

my $usage = "usage: $0 bwasw.sam.align_stats [longest_contig_only] [no_require_100pct_align]\n\n";


my $align_stats_file = $ARGV[0] or die $usage;
my $longest_contig_only_flag = $ARGV[1] || 0;
my $no_require_100pct_align = $ARGV[2] || 0;


=format

0       #
1       scaff_name
2       read_name
3       read_length
4       aligned_bases
5       matches
6       mismatches
7       indel_bkpts
8       sum_indel_lens
9       pct_mismatches
10      pct_indel_bkpts
11      pct_indel_lens


=cut


    ;


my @pct_mismatches;
my @pct_indels;

my $sum_length = 0;
my $sum_mismatches = 0;
my $sum_indels = 0;


my @perfect_aligns;
my @imperfect_aligns;

open (my $fh, $align_stats_file) or die $!;



my %core_acc_to_entries;


my $prev_read_name = "";

while (<$fh>) {
    chomp;
    my $line = $_;

    my @x = split(/\t/);
    if (scalar (@x) == 12 && $x[3] =~ /^\d+$/) {

        my $read_name = $x[2];
        if ($read_name eq $prev_read_name) { 
            next; # only one alignment per read
        }
        $prev_read_name = $read_name;
        


        my $core_read_name = $read_name;
        $core_read_name =~ s/_\d+$//;
        $core_read_name =~ s/\.\d+\.PbioCR$//;
        my $read_length = $x[3];
        
        push (@{$core_acc_to_entries{$core_read_name}}, { line => $line,
                                                          read_len => $read_length,
                                                      });
    }
}
close $fh;


foreach my $core_read_acc (keys %core_acc_to_entries) {

    my @alignments = @{$core_acc_to_entries{$core_read_acc}};

    if ($longest_contig_only_flag) {
        @alignments = sort {$a->{read_len}<=>$b->{read_len}} @alignments;
        my $longest_read = pop @alignments;
        @alignments = ($longest_read);
    }

    foreach my $alignment (@alignments) {


        my $line = $alignment->{line};
        my @x = split(/\t/, $line);
        

        my $seq_length = $x[3];
        my $aligned_bases = $x[4];
        my $mismatches = $x[6];
        my $indels = $x[8];
        
        my $pct_mism = $x[9];
        my $pct_indl = $x[11];
        
        
        push (@pct_mismatches, $pct_mism);
        push (@pct_indels, $pct_indl);
        
        $sum_length += $aligned_bases;
        $sum_mismatches += $mismatches;
        $sum_indels += $indels;
        
        
        if ($mismatches == 0 && $indels == 0 && ($no_require_100pct_align || $seq_length == $aligned_bases)) {
            push (@perfect_aligns, $line);
        }
        else {
            push (@imperfect_aligns, $line);
        }
        
    }
}


my $avg_pct_mismatch = $sum_mismatches / $sum_length * 100;
my $avg_pct_indel = $sum_indels / $sum_length * 100;

print "// base stats:\n";
print "Avg_pct_mismatch_all_bases: $avg_pct_mismatch\n";
print "Avg_pct_indel_all_bases: $avg_pct_indel\n";


print "\n// assembly stats:\n";
my $stat = Statistics::Descriptive::Sparse->new();
$stat->add_data(@pct_mismatches);
print "Mean pct_mismatch of assembly: " . $stat->mean() . "\n";

$stat->clear();
$stat->add_data(@pct_indels);
print "Mean pct_indel of assembly: " . $stat->mean() . "\n";
print "\n\n";



################################################
#   write some files for downstream analysis
################################################


# write files for interrogation using R
open (my $ofh, ">$align_stats_file.pct_mismatch.dat") or die $!;
print $ofh join("\n", @pct_mismatches) . "\n";
close $ofh;

open ($ofh, ">$align_stats_file.pct_indel.dat") or die $!;
print $ofh join("\n", @pct_indels) . "\n";
close $ofh;

open ($ofh, ">$align_stats_file.perfect_aligns.txt") or die $!;
print $ofh join("\n", @perfect_aligns) . "\n" if @perfect_aligns;
close $ofh;


open ($ofh, ">$align_stats_file.imperfect_aligns.txt") or die $!;
print $ofh join("\n", @imperfect_aligns) . "\n" if @imperfect_aligns;
close $ofh;

open ($ofh, ">$align_stats_file.all_aligns.txt") or die $!;
print $ofh join("\n", @perfect_aligns, @imperfect_aligns) . "\n" if (@perfect_aligns || @imperfect_aligns);
close $ofh;



exit(0);

