#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 file.intron_analysis ...\n\n";

my @files = @ARGV;
unless (@files) {
    die "$usage\n\n";
}


my @types = qw(
               NOVEL
               REF
               NOVELCOMBO
               REFCOMBO
               SUBrefCombo);

print "#analysis\t" . join("\t", @types) . "\n";

foreach my $intron_analysis_output (@files) {
        
    my %T = map { + $_ => 1 } @types;
    
    my %type_to_feature;
    
    open (my $fh, $intron_analysis_output) or die "Error, cannot open file $intron_analysis_output";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my ($type, $intron_set, @rest) = split(/\t/);
        
        unless ($T{$type}) {
            print STDERR "Error, do not understand type: $type\n$_\n";
            next;
        }
        
        $type_to_feature{$type}->{$intron_set} = 1;
    }
    close $fh;
    
    
    my @vals;
    foreach my $type (@types) {
        my $count = scalar (keys %{$type_to_feature{$type}});
        push (@vals, $count);
    }
    print $intron_analysis_output . "\t" . join("\t", @vals) . "\n";
    
}

exit(0);


