#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 transcripts.fasta jaccard_clips.wig\n\n";

my $transcript_fa = $ARGV[0] or die $usage;
my $jaccard_clips = $ARGV[1] or die $usage;



main: {

	my %trans_acc_to_clips = &get_clip_pts($jaccard_clips);
	
	my $fasta_reader = new Fasta_reader($transcript_fa);
	
	while (my $seq_obj = $fasta_reader->next()) {

		my $acc = $seq_obj->get_accession();
		my $header = $seq_obj->get_header();
        my $seq = $seq_obj->get_sequence();

		my @parts = split(/;/, $acc);
		my $kmer_cov = pop @parts;

        my ($header_begin, $rest_header) = split(/\s+/, $header, 2);
                
		if (my $clips_aref = $trans_acc_to_clips{$acc}) {
			
			##
			my $start = 1;
			while (@$clips_aref) {
				my $clip = shift @$clips_aref;
				my $length = $clip - $start + 1;
				my $subseq = substr($seq, $start-1, $length);
				print ">$acc.$start-$clip;$kmer_cov $rest_header\n$subseq\n";
				$start = $clip + 1;
			}
			if ($start < length($seq) + 25) {
                # require subseq to be at least 25 bases long (one kmer length)
                my $subseq = substr($seq, $start-1, length($seq)-$start + 1);
                print ">$acc.$start-" . length($seq) . ";$kmer_cov $rest_header\n$subseq\n"; # coverage value needs to be the last piece of the accession for use by PASA
			}
		}
		else {
            # no clippint
			print ">$header\n$seq\n";
		}
	}
    

	exit(0);

}


####
sub get_clip_pts {
	my ($jaccard_clips) = @_;

	my %trans_to_clips;

	my $trans_acc = "";

	open (my $fh, $jaccard_clips) or die "Error, cannot open file $jaccard_clips";
	while (<$fh>) {
		chomp;
		if (/^variableStep chrom=(\S+)/) {
			$trans_acc = $1;
		}
		elsif (/^(\d+)/) {
			my ($coord, $val) = split(/\t/);
            if ($val) {
                push (@{$trans_to_clips{$trans_acc}}, $coord);
            }
		}
	}
	close $fh;

	return(%trans_to_clips);
}


			
