#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../../PerlLib");
use Nuc_translator;

my $usage = "usage: $0 blat.psl.nameSorted.sam  reads.tab.nameSorted\n\n";


my $sam = $ARGV[0] or die $usage;
my $seqs = $ARGV[1] or die $usage;


main: {

	open (my $sam_fh, "$sam") or die "Error, cannot open file $sam";
	open (my $seqs_fh, "$seqs") or die "Error, cannot open file $seqs";
	
	my $seq_line = <$seqs_fh>;
	chomp $seq_line;
	my ($seq_acc, $seq, $qual) = split(/\t/, $seq_line);
	$seq_acc =~ s/\s//g; # rid any ws from acc name

    unless ($qual) {
        $qual = 'B' x length($seq);
    }

	while (my $sam_line = <$sam_fh>) {
		
		my @x = split(/\t/, $sam_line);
		
		my $acc = $x[0];
		my $flag = $x[1];

		my $aligned_orient = ($flag & 0x0010) ? '-' : '+';
		

		while ($seq_acc lt $acc) {
			$seq_line = <$seqs_fh>;
			chomp $seq_line;
			($seq_acc, $seq, $qual) = split(/\t/, $seq_line);
			$seq_acc =~ s/\s//g; # no ws in acc name
            unless (defined $qual) {
				$qual = 'B' x length($seq); #$qual = "*";
			}
		}
		
		if ($acc eq $seq_acc) {
	
			if ($aligned_orient eq '-') {
				my $revseq = &reverse_complement($seq);
				my @q = split(//, $qual);
				my $revqual = join("", reverse(@q));
				$x[9] = $revseq;
				$x[10] = $revqual;
			}
			else {
				$x[9] = $seq;
				$x[10] = $qual;
			}
		}
		else {
			die "Error,\n[$acc]\nnot encountered in file: $seqs,\ncurrently cursor is at seq:\n[$seq_acc]\n";
		}
		
		## set read quality score to a high value so that Scripture will use it.
		$x[4] = 255;
		$x[5] =~ s/H/S/gi; # convert hard to soft clips; samtools doesn't like H's
		
		foreach my $val (@x) {
			unless (defined $val) {
				$val = "*";
			}
		}
		
		print join("\t", @x);
	}

	exit(0);
}


