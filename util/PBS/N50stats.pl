#!/usr/bin/env perl

=pod

=head1 NAME

=head1 USAGE

	-in input file in FASTA/Q or posmap length file (e.g. posmap.scflen)
	-genome genome size in bp for estimating N lengths and indexes
	-single FASTA/Q has sequence in a single line (faster)
	-overwrite => Force overwrite
	-reads	=> Force processing as read data (no N50 statistics)
	-noreads => Force as not being read data. Good for cDNA assemblies with short contigs

=head1 AUTHORS

 Alexie Papanicolaou 1
	
	Ecosystem Sciences, CSIRO, Black Mountain Labs, Clunies Ross Str, Canberra, Australia
	alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

None known so far.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Statistics::Descriptive;
use Bio::SeqIO;
$|=1;

my (@infiles,$user_genome_size,$is_fasta,$is_fastq,$is_single,$overwrite,$is_reads,$isnot_reads);
GetOptions(
	'in=s{,}'    => \@infiles,
	'single' =>\$is_single,
	'genome:s' => \$user_genome_size,
	'overwrite' => \$overwrite,
	'reads'	=>\$is_reads,
	'noreads'	=>\$isnot_reads,
);
if (!@infiles){
	@infiles = @ARGV;
}
pod2usage "No input files!\n" if !@infiles;
die "Cannot ask for both reads and noreads options at the same time!\n" if $is_reads && $isnot_reads;
 
if ($is_reads && !$isnot_reads){
	print "Processing all data as reads\n";
}
if ($user_genome_size && $user_genome_size=~/\D/){
	if ($user_genome_size=~/^(\d+)kb$/i){
		$user_genome_size=int($1.'000');
	}
	elsif ($user_genome_size=~/^(\d+)mb$/i){
		$user_genome_size=int($1.'000000');
	}
	elsif ($user_genome_size=~/^(\d+)gb$/i){
		$user_genome_size=int($1.'000000000');
	}
	print "Genome set to ".&thousands($user_genome_size)." b.p.\n";
}

foreach my $infile (@infiles){
	unless ($infile && -s $infile){warn("I need a posmap length file, e.g. .posmap.scflen for scaffolds\n");pod2usage;}
	my $outfile=$infile.'.n50';
	$outfile.='g' if ($user_genome_size);
	warn ("Outfile $outfile already exists\n") if -s $outfile && !$overwrite;
	next  if -s $outfile && !$overwrite;
	my $total=int(0);
	my $gaps = int(0);
	my $seq_ref;
	my @head=`head $infile`;
	foreach (@head){
		if ($_=~/^>\S/){
			$is_fasta=1;
			print "FASTA file found!\n";
			last;
		}elsif($_=~/^@\S/){
			$is_fastq=1;
			print "FASTQ file found!\n";
			last;
		}
	}

	print "Parsing file $infile...\n";
	if ($is_fasta){
		($total,$gaps,$seq_ref) = &process_fasta($infile);
	}
	elsif($is_fastq){
		($total,$gaps,$seq_ref) = &process_fastq($infile);
	}
	else {
		($total,$gaps,$seq_ref) = &process_csv($infile);
	}
	print "Preparing stats...\n";
	my ($mean,$n50,$n10,$n25,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum,$genome_size) = &process_stats($seq_ref,$total);

	if ($mean){
		open (OUT,">".$outfile);
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data($seq_ref);
		#my $skew='';sprintf("%.2f",$stat->skewness());
		my $mean = sprintf("%.2f",$mean);
		my $median = $stat->median();
		my $var = sprintf("%.2f",$stat->variance());
		my $sd = sprintf("%.2f",$stat->standard_deviation());
		if (!$scaffolds || $scaffolds == 0){
			$scaffolds=$sequence_number;
			$scaffolds_size=$smallest;
		}
		print OUT "File: $infile\n";
		print OUT "TOTAL: ".&thousands($total)." bp in ".&thousands($sequence_number)." sequences\n";
		print OUT "\tof which ".&thousands($gaps)." are Ns/gaps.\n";
		print OUT "Mean: ".&thousands($mean)."\nStdev: ".&thousands($sd)."\n";
		print OUT "Median: ".&thousands($median)."\n";
		print OUT "Smallest: ".&thousands($smallest)."\nLargest: ".&thousands($largest)."\n";
		if (($mean >=1000 && !$is_reads) || $isnot_reads){
			print OUT "N10 length: ".&thousands($n10_length)."\nN10 Number: ".&thousands($n10)."\n";
			print OUT "N25 length: ".&thousands($n25_length)."\nN25 Number: ".&thousands($n25)."\n";
			print OUT "N50 length: ".&thousands($n50_length)."\nN50 Number: ".&thousands($n50)."\n";
			print OUT "Assuming a genome size of "
				.&thousands($user_genome_size)
				." then the top ".&thousands($scaffolds)
				." account for it (min "
				.&thousands($scaffolds_size)
				." bp)\n" if $user_genome_size;
		}else{
			print OUT "Reads found! Read coverage estimated to ".sprintf("%.2f",$total/$user_genome_size)."x using user provided genome size of ".&thousands($user_genome_size)."\n" if $user_genome_size;
		}
		close (OUT);
		print "Done, see $outfile\n";
		system("cat $outfile");
	}else {
		open (OUT,">".$outfile);
		print OUT "File: $infile\n";
		print OUT "TOTAL: $total bp in $sequence_number sequences\n";
		close (OUT);
		warn "Non fatal warning: Something went wrong in estimating the statistics. Maybe the provided genome length is much larger than sequence length or maybe less than 3 sequences provided?\n";
	}
}
########################################################################
sub process_fasta(){
	print "Processing as FASTA\n";
	my $infile=shift;
	my @array ;
	my $total=int(0);
	my $gaps=int(0);
	my $counter = int(0);

	if ($is_single){
		open (IN,$infile)||die;
		while (my $seq_id=<IN>) {
	                my $seq=<IN>;
	                my $length=length($seq)-1; # newline
			$counter+=length($seq_id)+$length+1;
	                next unless $length;
			$gaps+=($seq=~tr/[N\-]//);
	                push(@array,$length);
        	        $total+=$length;
		}
		close IN;
	}else{
		my $filein = new Bio::SeqIO(-file=>$infile , -format=>'fasta');
		while (my $seq_obj=$filein->next_seq()) {
			$counter+=length($seq_obj->seq().$seq_obj->description().' '.$seq_obj->id()) if $seq_obj->seq();
			my $length=$seq_obj->length();
			my $seq=$seq_obj->seq();
			$gaps+=($seq=~tr/[N\-]//);
			next unless $length;
			push(@array,$length);
	        	$total+=$length;
	        }
	}
	print "\n";
	die "No data found or wrong format\n" unless $total;
	return ($total,$gaps,\@array);
}
sub process_fastq(){
	print "Processing as FASTQ\n";
	my $infile=shift;
	my @array ;
	my $total=int(0);
	my $gaps=int(0);
	my $counter = int(0);
	open (IN,$infile);
	while (my $seq_id=<IN>) {
		my $seq=<IN>;
		my $scrap=<IN>.<IN>;
		my $length=length($seq)-1; #newline
		$counter+=length($seq_id)+$length+1;
		next unless $length;
		$gaps+=($seq=~tr/[N\-]//);
                push(@array,$length);
                $total+=$length;
       }
	close IN;
	print "\n";
	die "No data found or wrong format\n" unless $total;
	return ($total,$gaps,\@array);
}
sub process_csv(){
	print "Processing as CSV\n";
	my $infile=shift;
	my @array ;
	my $total=int(0);
	my $gaps='N/A';
	my $counter = int(0);
	open (IN,$infile)||die ("Cannot open $infile\n");
	while (my $ln=<IN>){
		$counter+=length($ln);
		$ln=~/(\d+)$/;
		next unless $1;
		my $length= $1;
		push(@array,$1);
		$total+=$length;
	}
	close IN;
	die "No data found or wrong format\n" unless $total;
	return ($total,$gaps,\@array);
}

sub process_stats(){
	my $sequences_ref = shift;
	my $total = shift;
	my $genome_size = int(0);
	my $mean = $total / scalar(@$sequences_ref);
	my ($n50,$n10,$n25,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum);
	print "Sorting...";
	my @sequences=sort{$b<=>$a} @$sequences_ref;
	$smallest=$sequences[-1];
	$largest=$sequences[0];
	print " done!\n";
	$|=0;
	if (($mean < 1000 && !$isnot_reads) || $is_reads ){
		print "Reads detected. Ignoring N* calculations.\n";
		$genome_size=$total;
		$sequence_number = scalar(@$sequences_ref);
		return ($mean,$n50,$n10,$n25,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum) if $mean <1000;
	}
	elsif (!$user_genome_size){
		print "Setting genome size for N* calculations to total consensus $total\n";
		$genome_size=$total;
	}elsif($user_genome_size){
		print "Setting genome size for N* calculations to user defined $user_genome_size\n";
		$genome_size = $user_genome_size;
	}
	
	foreach my $sequence_length ( @sequences){
		$sum+=$sequence_length;
		$sequence_number++;
		if($sum >= $genome_size*0.1 && !$n10){
			$n10=$sequence_number;
			$n10_length=$sequence_length;
		}
		elsif($sum >= $genome_size*0.25 && !$n25){
			$n25=$sequence_number;
			$n25_length=$sequence_length;
		}
		elsif($sum >= $genome_size*0.5 && !$n50){
			$n50 = $sequence_number;
			$n50_length=$sequence_length;
		}elsif ($sum >= $genome_size && !$scaffolds){
			$scaffolds = $sequence_number;
			$scaffolds_size = $sequence_length;
		}

	}
	print "Processed $sequence_number sequences\n";
	return ($mean,$n50,$n10,$n25,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum,$genome_size);
}

sub thousands($){
	my $val = shift;
	return int(0) if !$val;
	$val = sprintf("%.0f", $val);
	1 while $val =~ s/(.*\d)(\d\d\d)/$1,$2/;
	return $val;
}
