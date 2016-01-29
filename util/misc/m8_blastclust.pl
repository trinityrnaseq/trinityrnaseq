#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Fasta_reader;
use SingleLinkageClusterer;
use File::Basename;
use Process_cmd;

my $usage = <<__EOUSAGE__;

########################################################################################
#
#  -I <string>        input fasta file
#  --prot or --nuc    protein or nucleotide
#
#  -L <int>           min percent length cutoff (default: 50)
#  -R                 the above percent length cutoff must be reciprocal (default: off)
#  -P <int>           min percent identity (default: 75)
#  
#  -E <float>         max E-value (default: 1e-10)
#  -H <int>           number of hits per entry (default: 100)
#
#  --megablast        use megablast (only for nuc mode)
#
#  --use_m8 <string>    use an existing blast m8 file  
#
#########################################################################################

__EOUSAGE__

	;



main: {

	my $fasta_file;
	my $prot_mode = 0;
	my $nuc_mode = 0;
	my $min_percent_length = 50;
	my $reciprocal_length_flag = 0;
	my $min_percent_identity = 75;
	my $Evalue = 1e-10;
	my $num_hits = 100;
	my $megablast_flag = 0;
	my $blast_m8_file = "";
    
	my $help_flag;
	
	
	&GetOptions ( 'h' => \$help_flag,
				  'I=s' => \$fasta_file,
				  'prot' => \$prot_mode,
				  'nuc' => \$nuc_mode,
				  'L=i' => \$min_percent_length,
				  'R' => \$reciprocal_length_flag,
				  'P=i' => \$min_percent_identity,
				  'E=f' => \$Evalue,
				  'H=i' => \$num_hits,
				  'megablast' => \$megablast_flag,
				  'use_m8=s' => \$blast_m8_file, 
                  
				  );
	
	
	
	
	if ($help_flag) {
		die $usage;
	}
	
	
	unless ($fasta_file && ($prot_mode || $nuc_mode) ) {
		die $usage;
	}
	
		
	
	my %seq_lengths = &parse_seq_lengths($fasta_file);
	
	$num_hits++; # include self hit
	
	my $blast_outfile = $blast_m8_file;
    unless ($blast_outfile) {

        $blast_outfile = join (".", "blast", $Evalue, $num_hits, "m8");
    }

    
	unless (-s $blast_outfile) {
	
		## prep fasta file for blast search and run blast
		if ($prot_mode) {
			unless (-s "$fasta_file.pin") {
				my $cmd = "formatdb -i $fasta_file -p T";
				&process_cmd($cmd);
			}
			
			my $cmd = "blastall -p blastp -d $fasta_file -i $fasta_file -m 8 -e $Evalue -v $num_hits -b $num_hits -F \'m S\' > $blast_outfile";
			&process_cmd($cmd);
			
		}
		else {
			#nuc mode
			unless (-s "$fasta_file.nin") {
				my $cmd = "formatdb -i $fasta_file -p F";
				&process_cmd($cmd);
			}
			
			my $prog = ($megablast_flag) ? "megablast" : "blastall -p blastn";

			my $cmd = "$prog -d $fasta_file -i $fasta_file -m 8 -e $Evalue -v $num_hits -b $num_hits -F \'m D\' > $blast_outfile";
			&process_cmd($cmd);
		}
	}
	
	## assign length and percent length coverage info 
	my $length_outfile = "$blast_outfile.wLens";
	unless (-s $length_outfile) {
		&compute_length_coverage($blast_outfile, $length_outfile, \%seq_lengths);
	}

	&cluster_hits($length_outfile, $min_percent_length, $min_percent_identity, $reciprocal_length_flag);

	exit(0);

}

####
sub cluster_hits {
	my ($length_outfile, $min_percent_length, $min_percent_identity, $reciprocal_length_flag) = @_;

	my @pairs;

	open (my $fh, $length_outfile) or die "Error, cannot read file $length_outfile";
	while (<$fh>) {
		chomp;
		my @x = split(/\t/);
		my $accA = $x[0];
		my $accB = $x[1];
		my $per_ID = $x[2];
		my $length_covA = $x[13];
		my $length_covB = $x[15];
		
		if ($accA eq $accB) { next; } # no self matches examined.

		if ($per_ID < $min_percent_identity) { next; }

		my $covA_requirements = ($length_covA >= $min_percent_length) ? 1:0;
		my $covB_requirements = ($length_covB >= $min_percent_length) ? 1:0;

		if ($reciprocal_length_flag) { 
			
			if ($covA_requirements && $covB_requirements) {
				push (@pairs, [$accA, $accB]);
			}
		}
		else {
			if ($covA_requirements || $covB_requirements) {
				push (@pairs, [$accA, $accB]);
			}
		}
	}
	close $fh;


	my $clusters_outfile = "$length_outfile.L$min_percent_length.P$min_percent_identity.R$reciprocal_length_flag.clusters";
	open (my $ofh, ">$clusters_outfile") or die "Error, cannot write to file $clusters_outfile";
	
	my @clusters = &SingleLinkageClusterer::build_clusters(@pairs);

	my %size_counter;
	foreach my $cluster (@clusters) {
		print $ofh join("\t", @$cluster) . "\n";
		
		my $num_eles = scalar(@$cluster);
		$size_counter{$num_eles}++;

	}
	close $ofh;
	
	## provide summary report for dist
	print "#cluster_size\tcount\n";
	foreach my $cluster_size (sort {$a<=>$b} keys %size_counter) {
		my $count = $size_counter{$cluster_size};
		print "$cluster_size\t$count\n";
	}
	
	return;
	
}




####
sub compute_length_coverage {
	my ($blast_results_file, $length_outfile, $seq_lengths_href) = @_;

	open (my $ofh, ">$length_outfile") or die "Error, cannot write to $length_outfile";
	open (my $fh, "$blast_results_file") or die "Error, cannot read file $blast_results_file";
	
	while (<$fh>) {
		chomp;
		my @x = split(/\t/);
		my ($accA, $end5_A, $end3_A) = ($x[0], $x[6], $x[7]);
		my ($accB, $end5_B, $end3_B) = ($x[1], $x[8], $x[9]);

		my $seq_lenA = $seq_lengths_href->{$accA} or die "Error, no sequence length for acc: $accA";
		my $seq_lenB = $seq_lengths_href->{$accB} or die "Error, no sequence length for acc: $accB";

		my $hit_lenA = abs($end3_A - $end5_A) + 1;
		my $hit_lenB = abs($end3_B - $end5_B) + 1;

		my $percent_lenA = sprintf("%.1f", $hit_lenA / $seq_lenA * 100);
		my $percent_lenB = sprintf("%.1f", $hit_lenB / $seq_lenB * 100);

		print $ofh join("\t", @x, $seq_lenA, $percent_lenA, $seq_lenB, $percent_lenB) . "\n";

	}
	
	close $ofh;

	return;
}





####
sub parse_seq_lengths {	
	my ($fasta_file) = @_;

	my %lengths;

	my $seq_lengths_file = basename($fasta_file) . ".seqLen";

	if (-e $seq_lengths_file) {
		open (my $fh, $seq_lengths_file) or die "Error, cannot open file $seq_lengths_file";
		while (<$fh>) {
			chomp;
			my ($len, $acc, @rest) = split(/\s+/);
			$lengths{$acc} = $len;
		}
		close $fh;

		return(%lengths);
	}
	else {
			
		my $fasta_reader = new Fasta_reader($fasta_file);
		
		open (my $ofh, ">$seq_lengths_file") or die "Error, cannot write to file $seq_lengths_file";
		
		while (my $seq_obj = $fasta_reader->next()) {
			my $acc = $seq_obj->get_accession();
			my $sequence = $seq_obj->get_sequence();
			
			my $len = length($sequence);
			$lengths{$acc} = $len;

			print $ofh "$len\t$acc\n";
		

		}
		close $ofh;

		return(%lengths);
	}

}
