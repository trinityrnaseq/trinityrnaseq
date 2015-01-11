#!/usr/bin/perl -w

# written by
#
# Bob Freeman, Ph.D.
# Acorn Worm Informatics, Kirschner lab
# Dept of Systems Biology, Alpert 524
# Harvard Medical School
# 200 Longwood Avenue
# Boston, MA  02115
# 617/432.2294, vox
#
# bob_freeman@hms.harvard.edu

=from_Bob

I wrote a small script in Perl that gives a number of basic stats and can even generate a histogram of lengths for you. Script is attached. I typically use this for assessing the size / quality of my data. 

For example, after trimming a set of reads (prior to assembly), I used the command:

    fastq_stats.pl -i trimmed_reads.fastq -f "30,40,50,60,70,80,90,100"

to assess resulting data, and my output is:

    Count27717551
    Sum2347820942
    Mean84.7052087

    Min36
    Max101
    Median101

    Q036
    Q154
    Q2101
    Q3101
    Q4101

    36(min)
    <=  30.0:0
    <=  40.0:183537
    <=  50.0:283964
    <=  60.0:8209709
    <=  70.0:339854
    <=  80.0:458850
    <=  90.0:460218
    <= 100.0:1271776

I believe I've modified the file to identify fasta, fastq, and protein fasta ( *.fasta, *.fa, *.mfasta, *.pro, and *.pep). Output can be in table-format also ( -t ). And use the -h switch to print out the (scarce) help info. It does use Bioperl and the Perl modules Statistics::Descriptive.

Hope you find it useful!

Bob

=cut


use strict;
use Bio::SeqIO;
use Carp;
use Data::Dumper;
use English	qw ( -no_match_vars );
use Statistics::Descriptive;


my $usage 							= qq (
	fastq_stats.pl -i infile -f partition|bins -t

	reads in fastq or fasta data and spits out some descripting stats
	-i	input file
	-f	frequency distribution, either # partitions or a quoted list of
		comma-separated bins
	-t  print results in table (tab-delimited) format

);

# opens fastq or fasta data
# reads in each sequence, grabbing data about each
# then spits out output data

$OUTPUT_AUTOFLUSH = 1;
our @seqlen_data;
our @hist_bins;
our ($infile, $freqdist, $fdref);
our $table = 0;

parse_args();
read_seq_data2();
if (!$table) {
    output_stats();
} else {
    output_stats_table();
}

exit;

#
#
# END OF PROGRAM
#

sub parse_args {
	#
	while (scalar(@ARGV)) {
		my $arg						= shift(@ARGV);
		if	  ($arg eq '-h') 		{ die 	$usage;             } 
		elsif ($arg eq '-i') 		{ $infile 	= shift(@ARGV); }
		elsif ($arg eq '-f') 		{ $freqdist = shift(@ARGV); }
		elsif ($arg eq '-t') 		{ $table    = 1; }
		else 						{ die 		"unknown argument '$arg'"; }
	}

	# check to ensure we have the required parameters
	if (!defined $infile) {
		print "Error from undefined input arguments!\n";
		die $usage;
	}	
	
	# parse frequency distribution stuff
	if (defined $freqdist) {
		$freqdist =~ s/ //g;
		@hist_bins = split (/,/, $freqdist);
		$freqdist = 1;
	} else {
		$freqdist = 0;	
	}
}

sub read_seq_data {

	my ($file_suffix) = $infile =~ m/.+\.(\S+)$/;
	# create input SeqIO object
	my $seq_in = Bio::SeqIO->new(	'-file'   => "<$infile",
                             		'-format' => $file_suffix );
                             		
    while (my $inseq = $seq_in->next_seq()) {
    	push @seqlen_data, $inseq->length();
    
    } #end of while
}

sub read_seq_data2 {
	
	if ($infile =~ m/(\.fasta$)|(\.fa$)|(\.mfasta$)|(\.mfa$)|(\.pro$)|(\.pep$)/) {
		read_fasta();
	} elsif ($infile =~ m/(\.fastq$)|(\.fq$)/) {
		read_fastq();	
	} else {
		die "Error: Unknown file format for input file $infile!\n";
	}
}

sub read_fasta {
	# read in fasta data using BioPerl methods. These are pretty quick.
	# create input SeqIO object
	my $seq_in = Bio::SeqIO->new(	'-file'   => "<$infile",
                             		'-format' => 'fasta' );
                             		
    while (my $inseq = $seq_in->next_seq()) {
    	push @seqlen_data, $inseq->length();
    
    } #end of while
}

sub read_fastq {
	# read in fastq data in old-fashioned method, as it's faster than using BioPerl
	open (my $INFILE, "<$infile")	|| die "Error opening file $infile for reading: $!\n";
	while ( my $line = <$INFILE> ) {
		#chomp $line;
		#next if $line =~ m/^\s|#/;
		if ($line =~ m/^@/) {
			# skip to next line and push length
			$line = <$INFILE>;
			chomp $line;
			push @seqlen_data, length($line);
		}
	}
	close $INFILE;
}

sub output_stats {
	
	print "\nSequence stats:\n";
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@seqlen_data);
	print "Count\t", $stat->count(), "\n";
	print "Sum\t", $stat->sum(), "\n";
	print "Mean\t", $stat->mean(), "\n\n";
	
	print "Min\t", $stat->min(), "\n";
	print "Max\t", $stat->max(), "\n";
	print "Median\t", $stat->median(), "\n\n";
	print "Q0\t", $stat->quantile(0), "\n";
	print "Q1\t", $stat->quantile(1), "\n";
	print "Q2\t", $stat->quantile(2), "\n";
	print "Q3\t", $stat->quantile(3), "\n";
	print "Q4\t", $stat->quantile(4), "\n";
	
	# do frequency distribution output if indicated
	if ($freqdist) {
		# call proper method
		if (scalar(@hist_bins) == 1) {
			$fdref = $stat->frequency_distribution_ref($hist_bins[0]);
		} else {
			$fdref = $stat->frequency_distribution_ref(\@hist_bins);
		}
		# now print it out
		print "\n", $stat->min(), "(min)\n";		
		for (sort {$a <=> $b} keys %$fdref) {
			#sprintf "<= %.1f: %8u \n", $_, "<= $_, \tcount = $fdref->{$_}\n";
			printf "<= %5.1f:\t%8u \n", $_, $fdref->{$_};
		}
	}
}

sub output_stats_table {
	#
	# do calculations
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@seqlen_data);
	if ($freqdist) {
		# call proper method
		if (scalar(@hist_bins) == 1) {
			$fdref = $stat->frequency_distribution_ref($hist_bins[0]);
		} else {
			$fdref = $stat->frequency_distribution_ref(\@hist_bins);
		}
	}

	# print header
	print "#Sequence stats:\n";
	print join("\t", "#","Count","Sum","Mean","Min","Max","Median","Q0","Q1","Q2","Q3","Q4");
    if ($freqdist) {	
		for (sort {$a <=> $b} keys %$fdref) {
			#sprintf "<= %.1f: %8u \n", $_, "<= $_, \tcount = $fdref->{$_}\n";
			print "\t<=", $_;
		}
    }
    print "\n";	
    
    #print results
    print $infile, "\t";
	print $stat->count(), "\t";
	print $stat->sum(), "\t";
	print $stat->mean(), "\t";
	print $stat->min(), "\t";
	print $stat->max(), "\t";
	print $stat->median(), "\t";
	print $stat->quantile(0), "\t";
	print $stat->quantile(1), "\t";
	print $stat->quantile(2), "\t";
	print $stat->quantile(3), "\t";
	print $stat->quantile(4);
	if ($freqdist) {
		# now print it out
		for (sort {$a <=> $b} keys %$fdref) {
			#sprintf "<= %.1f: %8u \n", $_, "<= $_, \tcount = $fdref->{$_}\n";
			print "\t", $fdref->{$_};
		}
	}
	print "\n";
}
