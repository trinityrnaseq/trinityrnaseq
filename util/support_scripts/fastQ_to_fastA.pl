#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Nuc_translator;
use IO::Uncompress::Gunzip;

use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

##########################################################
#
#  -I <string>     input.fq  or "input1.fq,input2.fq,input3.fq"
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

&GetOptions( 'I=s'          => \$inputFile,
             'ignore_dirty' => \$ignore_dirty,
             'a=i'          => \$append_num,
             'rev'          => \$revcomp_flag,
	);


unless ($inputFile) {
	die $usage;
}

main: {
    my @files = split(/,/, $inputFile);
    foreach my $file (@files) {
        $file =~ s/\s//g;
        if ($file =~ /\w/) {
            &fastQ_to_fastA($file);
        }
    }
    exit(0);
}


sub fastQ_to_fastA {
    my ($file) = @_;

    my $fh = new IO::Uncompress::Gunzip($file) or die "Error, cannot open file $file";

    my $counter = 0;
    my $num_clean = 0;
    my $num_dirty = 0;

    while (my $line = <$fh>) {
        $line =~ s/\cM//g; # remove any cntrl-M characters (sometimes derived from MS-windows text files)

        if ($line =~ /^\@/) {
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

            if ($header && $seq && $qual_header && $qual_line =~ /\S/ &&
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
                if ($revcomp_flag) {
                    $seq = &reverse_complement($seq);
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
    }

    return;
}
