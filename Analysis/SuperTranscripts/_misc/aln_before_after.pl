#!/usr/bin/env perl

use strict;
use warnings;

use lib($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "\n\n\tusage: $0 trin_before.fa  trin_after.fa\n\n";

my $trin_before_file = $ARGV[0] or die $usage;
my $trin_after_file = $ARGV[1] or die $usage;

main: {

    my %trin_seqs_before = &parse_trin_fa($trin_before_file);

    my %trin_seqs_after = &parse_trin_fa($trin_after_file);


    foreach my $acc (keys %trin_seqs_after) {
        my $seq_before = $trin_seqs_before{$acc} or die "Error, no before transcript for [$acc]";
        my $seq_after = $trin_seqs_after{$acc} or die "Error, no after transcript for [$acc]";
        
        if ($seq_before eq $seq_after) {
            print "* $acc identical\n";
        }
        else {
            my $filename = $acc;
            $filename =~ s/\W/_/g;

            open(my $ofh, ">$filename") or die $!;
            print $ofh ">before\n$seq_before\n"
                . ">after\n$seq_after\n";
            close $ofh;

            my $ret = system("clustalw2 $filename");
            if ($ret) {
                die $!;
            }
        }
    }

    exit(0);

}

####
sub parse_trin_fa {
    my ($file) = @_;

    my $fasta_reader = new Fasta_reader($file);

    my %seqs = $fasta_reader->retrieve_all_seqs_hash();

    return(%seqs);
}
