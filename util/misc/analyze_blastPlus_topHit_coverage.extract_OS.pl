#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

#########################################################################################
#
# Required:
#
#  --blast_outfmt6_w_pct_hit_length <string>     blast.outfmt6.w_pct_hit_length  
#                                                (results from running analyze_blastPlus_tophat_coverage.pl)
#  
#  Optional:
#
#  --min_pct_hit_length <int>                    minimum percent hit length to be included in analysis. (default: 20)
#
#  --min_pct_species_report <float>              minimum percent of total species content to 
#                                                be reported in output (default: 1.0)
#
###########################################################################################


__EOUSAGE__

    ;

my $file;
my $min_pct_hit_len = 200;
my $min_pct_species_report = 1.0;

my $help_flag;

&GetOptions( 'blast_outfmt6_w_pct_hit_length=s' => \$file,
             'min_pct_hit_length=i' => \$min_pct_hit_len,
             'min_pct_species_report=f' => \$min_pct_species_report,
             
             'help|h' => \$help_flag,
             );


if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, didn't parse parameters: @ARGV";
}

unless ($file) {
    die $usage;
}


my $total = 0;
my %OS_counter;

my $prev_acc = "";

open (my $fh, $file) or die $!;
while (<$fh>) {
    if (/^\#/) { next; }
    
    chomp;
    my @x = split(/\t/);

    my $acc = $x[0];
    if ($acc eq $prev_acc) { next; }
                         

    my $pct_hit_len = $x[13];
    if ($pct_hit_len < $min_pct_hit_len) {
        next;
    }

    my $annot = $x[14];
    
    if ($annot =~ /OS=(\S+)/) {
        my $os = $1;
        $total++;
        $OS_counter{$os}++;
        
        $prev_acc = $acc;

    }
}
close $fh;

foreach my $os (sort {$OS_counter{$b} <=> $OS_counter{$a}} keys %OS_counter) {
    
    my $count = $OS_counter{$os};
    my $pct = sprintf("%.2f", $count/$total*100);

    print join("\t", $os, $count, "$pct%") . "\n" if $pct >= $min_pct_species_report;
}

exit(0);


