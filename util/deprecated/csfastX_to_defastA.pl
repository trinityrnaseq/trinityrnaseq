#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Nuc_translator;

use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

##########################################################
#
#  -I <string>     input.cfa/cfq
#
#  --ignoreDirty   ignores poorly formed entries
#
#  -a <int>        append "/num" to the accession name.
#
#  --rev           reverse complement nucleotide sequence.
#
###########################################################

_EOUSAGE_

	;

my $inputFile;
my $ignore_dirty = 0;
my $append_num;
my $revcomp_flag = 0;

&GetOptions( 'I=s' => \$inputFile,
			 'ignore_dirty' => \$ignore_dirty,
			 'a=i' => \$append_num,
			 'rev' => \$revcomp_flag,
	);


unless ($inputFile) {
	die $usage;
}

main: {

    my @files = split(/,/, $inputFile);
    
    foreach my $file (@files) {
        $file =~ s/\s//g;
        if ($file =~ /\w/) {
            &csfastX_to_defastA($file);
        }
    }
    

    exit(0);
}


sub csfastX_to_defastA(){
 my $infile = shift;
 open (my $fh, $infile) or die "Error, cannot open $infile";

 my $counter = 0;
 my $num_clean = 0;
 my $num_dirty = 0;

 my $line = <$fh>;

 while ($line) {
	if($line =~ /^#/){
		# skip over comments (present in some csfasta/csfastq files)
		$line = <$fh>;
		next;
	}

	elsif(substr($line,0,1) eq ">") {
		print STDERR "\r[$counter] [$num_clean clean] [$num_dirty dirty]       " if ($counter % 100000 == 0);

		# assume fasta format
		$counter++;

		my $header = $line;
		substr($header,0,1,''); # strip beginning ">"
		$line = <$fh>;
		my $seq = "";

		while($line && (substr($line,0,1) ne ">")){
			chomp $line;
			$seq .= $line;
			## find sequence
			$line = <$fh>;
		}

		my @header_parts = split(/\s+/, $header);
		$header = shift @header_parts;
		if (@header_parts && $header !~ m|/[12]$| && $header_parts[0] =~ /^([12])\:/) {
			$header .= "/$1";
		}
		
		if (defined $append_num) {
			$header .= "/$append_num";
		}
		
		if((substr($seq,0,1) =~ tr/ACGT//) > 0){
			# trim off primer base and first color (which depends on primer base)
			substr($seq,0,2,'');
		}
		# convert to double-encoded sequence
		$seq =~ tr/0123./ACGTN/;

		if ($revcomp_flag) {
			# reverse complement in color-space is just reverse
			$seq = reverse($seq);
		}

		print ">$header\n$seq\n";
		$num_clean++;

		# return to top of loop ($line should start with '>')
		next;
	}

	elsif(substr($line,0,1) eq "@") {
		# assume fastq format
		$counter++;
		
		# print STDERR "\r[$counter] [$num_clean clean] [$num_dirty dirty]       " if ($counter % 10000 == 0);

		my $header = $line;
		my $seq = <$fh>;
		my $qual_header = <$fh>;
		my $qual_line   = <$fh>;
				
		chomp $header;
		chomp $seq if $seq;
		chomp $qual_header if $qual_header;
		chomp $qual_line if $qual_line;
		
		if ($header && $seq && $qual_header && $qual_line && 
			$qual_header =~ /^\+/ && length($seq) == length($qual_line)) {
			
			# can do some more checks here if needed to be sure that the lines are formatted as expected.

			substr($header,0,1,''); # strip beginning "@"

			my @header_parts = split(/\s+/, $header);
			$header = shift @header_parts;
			if (@header_parts && $header !~ m|/[12]$| && $header_parts[0] =~ /^([12])\:/) {
			    $header .= "/$1";
			}
						

			if (defined $append_num) {
				$header .= "/$append_num";
			}
			
			# check for primer base at start
			if((substr($seq,0,1) =~ tr/ACGT//) > 0){
				# trim off primer base and first color (which depends on primer base)
				substr($seq,0,2,'');
			}

			# convert to double-encoded sequence
			$seq =~ tr/0123./ACGTN/;

			if ($revcomp_flag) {
				# reverse complement in color-space is just reverse
				$seq = reverse($seq);
			}
			
			print ">$header\n$seq\n";
			$num_clean++;
		}
		else {

			$num_dirty++;

			unless ($ignore_dirty) {
				die "Error, improperly formatted entry:\n\n$header\n$seq\n$qual_header\n$qual_line\n";
			}
			
		}
	}
	$line = <$fh>;
 }
}

