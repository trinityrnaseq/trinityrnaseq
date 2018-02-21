#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;


unless ($ARGV[0] && $ARGV[1]) {
    die "Usage: $0 \$fastaFile \$acc_listing\n";
}


my $fastaFile = $ARGV[0];
my $acc_listing = $ARGV[1];

open (ACC, "$acc_listing");
my %accs;
while (<ACC>) {
    chomp;
    my @x = split (/\s+/);
    foreach my $acc (@x) {
        if ($acc =~ /\w/) {
            $acc =~ s/\s//g;
            $acc =~ s/\W/_/g;
            $accs{$acc} = 1;
        }
    }
}
close ACC;



open (FASTA, "$fastaFile");

my %to_delete = %accs;

my $ok_flag = 1;
while (<FASTA>) {
    if (/^>(\S+)/) {
        my $acc = $1;
        if ($accs{$acc}) {
            # targeted for skipping
            $ok_flag = 0;
            delete $to_delete{$acc};
            print STDERR "-found and skipping $acc\n";
        }
        else {
            # not targeted for skipping
            $ok_flag = 1;
        }
    }
    if ($ok_flag) {
        print;
    }
}

if (%to_delete) {
    confess "Error, didn't observe and exclude the following accessions from the fasta file: " . Dumper(\%to_delete);
}

print STDERR "-done.\n\n";

exit(0);    
