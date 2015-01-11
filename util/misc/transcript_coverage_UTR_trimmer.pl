#!/usr/bin/env perl

use strict;
use warnings;

use threads;

use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);
use lib ("$FindBin::Bin/../../PerlLib");
use WigParser;
use Fasta_reader;
use Statistics::Descriptive;

my $usage = <<_EOUSAGE_;

########################################################################################################
#
#  Required:
#
#  --coord_sorted_SAM <string>      coordinate-sorted SAM file.
#  --transcripts <string>           transcripts in fasta file format
#
#  *If Strand-specific, specify:
#  --SS_lib_type <string>          library type:  if single: F or R,  if paired:  FR or RF
#  --trim_pct_median <int>         percent of median coverage value to set as upper threshold for end-trimming (default: 10)
#
#  --no_trim                       just provide the coverage info
#
########################################################################################################


_EOUSAGE_

	;


my $help_flag;

#my $partition_join_size;
my $SAM_file;
my $SS_lib_type = "";
my $transcripts_fasta_file = "";

my $trim_pct_median = 10;
my $NO_TRIM = 0;

&GetOptions ( 'h' => \$help_flag,

			  'coord_sorted_SAM=s' => \$SAM_file,
			  'SS_lib_type=s' => \$SS_lib_type,
			  'transcripts=s' => \$transcripts_fasta_file,
              'trim_pct_median=i' => \$trim_pct_median,
              'no_trim' => \$NO_TRIM,
              );


if ($help_flag) {
	die $usage;
}

unless (
	$SAM_file
    &&
    $transcripts_fasta_file
		) {
	die $usage;
}

if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
	die "Error, invalid --SS_lib_type, only F, R, FR, or RF are possible values";
}

my $UTIL_DIR = "$FindBin::Bin/";

main: {

	my $sam_file_to_process;

	if ($SS_lib_type) {
		
        $sam_file_to_process = "$SAM_file.+.sam";
        my $cmd = "$UTIL_DIR/SAM_strand_separator.pl $SAM_file $SS_lib_type";
        &process_cmd($cmd) unless (-s $sam_file_to_process);
		
	}
	else {
		$sam_file_to_process = $SAM_file;
	}
    
    ## define fragments
    my $cmd = "$UTIL_DIR/SAM_to_frag_coords.pl --sam $sam_file_to_process --min_insert_size 1 --max_insert_size 1000 "; ## writes file: $sam_file.frag_coords
    &process_cmd($cmd) unless (-s "$sam_file_to_process.frag_coords");
    
    ## define coverage
    $cmd = "$UTIL_DIR/fragment_coverage_writer.pl $sam_file_to_process.frag_coords > $sam_file_to_process.frag_coverage.wig";
    &process_cmd($cmd) unless (-s "$sam_file_to_process.frag_coverage.wig");
    
    my $wig_parser = new WigParser("$sam_file_to_process.frag_coverage.wig");
    
    my $fasta_reader = new Fasta_reader($transcripts_fasta_file);
    
    while (my $seq_obj = $fasta_reader->next()) {
        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();

        
        my $trim_5p_pos = "NA";
        my $trim_3p_pos = "NA";
        my $median = 0;
        my @cov = ();
        
        eval {
            # throws exception if there's no read coverage
            
            @cov = $wig_parser->get_wig_array($acc);
            #print ">$acc\n$sequence\n" . join(" ", @cov) . "\n";
            
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@cov);
            $median = $stat->median();
            #print "Median: $median\n";
            
            my $trim_thresh = int($trim_pct_median/100 * $median + 0.5);
            if ($trim_thresh == 0) {
                $trim_thresh = 1;
            }
            #print "Trim_thresh: $trim_thresh\n";
            

      
            
            unless ($NO_TRIM) {
                
                $trim_5p_pos = &trim_5p(\@cov, $trim_thresh);
                $trim_3p_pos = &trim_3p(\@cov, $trim_thresh);
            }
        };
        
        print join("\t", $acc, "$trim_5p_pos-$trim_3p_pos", $median, $sequence, join(" ", @cov) ) . "\n";
        
    }
    
    exit(0);
    
    
}


####
sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";
	
	my $ret = system($cmd);

	if ($ret) {
		die "Error, command $cmd died with ret $ret";
	}

	return;
}


####
sub trim_5p {
    my ($cov_aref, $threshold) = @_;

    for (my $i = 0; $i <= $#$cov_aref; $i++) {
        my $cov = $cov_aref->[$i];
        if ($cov >= $threshold) {
            return($i+1);
        }
    }

    return(-1); # bad sequence
    
}

####
sub trim_3p {
    my ($cov_aref, $threshold) = @_;
    
    for (my $i = $#$cov_aref; $i >=0; $i--) {
        my $cov = $cov_aref->[$i];
        if ($cov >= $threshold) {
            return($i+1);
        }
    }
    
    return(-1); # bad sequence
    
}
