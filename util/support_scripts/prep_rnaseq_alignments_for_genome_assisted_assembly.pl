#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use threads;

use File::Basename;
use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);
use lib("$FindBin::Bin/../../PerlLib");
use Thread_helper;


my $CPU = 2;

my $usage = <<_EOUSAGE_;

########################################################################################################
#
#  Required:
#
#  --coord_sorted_SAM <string>      coordinate-sorted SAM file.
#
#  -I  <int>                       maximum intron length  
#                                     (reads with longer intron lengths are ignored, and fragment reads 
#                                      farther apart on the genome are treated as unpaired))
#  -C  <int>                       min coverage for region boundary (default: 1)
#                             
#  *If Strand-specific, specify:
#  --SS_lib_type <string>          library type:  if single: F or R,  if paired:  FR or RF
#
#
#  Optional:
#
#  --min_reads_per_partition <int>      default: 10 
#  --parts_per_directory <int>          default: 100
#
#  --sort_buffer <string>               default: '10G'  amount of RAM to allocate to sorting.
#
#  --CPU <int>                        number of threads
#
########################################################################################################


_EOUSAGE_

	;

#  -J  <int>                       region join length (neighboring coverage bins within this range are merged into larger piles)



my $help_flag;

#my $partition_join_size;
my $max_intron_length;
my $SAM_file;
my $SS_lib_type = "";
my $min_coverage = 1;
my $min_reads_per_partition = 10;
my $parts_per_dir = 100;
my $sort_buffer = '10G';


&GetOptions ( 'h' => \$help_flag,

			  'coord_sorted_SAM=s' => \$SAM_file,
			  'SS_lib_type=s' => \$SS_lib_type,
			  
			  #'J=i' => \$partition_join_size,
			  'I=i' => \$max_intron_length,
              'C=i' => \$min_coverage,
              
              'min_reads_per_partition=i' => \$min_reads_per_partition,
              'parts_per_directory=i' => \$parts_per_dir,
              'CPU=i' => \$CPU,

              'sort_buffer=s' => \$sort_buffer,
              );


if ($help_flag) {
	die $usage;
}

unless (
	$SAM_file
        # && $partition_join_size
	&& $max_intron_length
		) {
	die $usage;
}

if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
	die "Error, invalid --SS_lib_type, only F, R, FR, or RF are possible values";
}

my $UTIL_DIR = "$FindBin::RealBin/";

main: {



    
	my @sam_info;

	if ($SS_lib_type) {
		my ($plus_strand_sam, $minus_strand_sam) = ("$SAM_file.+.sam", "$SAM_file.-.sam");
		if (-s $plus_strand_sam && $minus_strand_sam) {
			print STDERR "-strand partitioned SAM files already exist, so using them instead of re-creating them.\n";
		}
		else {
			my $cmd = "$UTIL_DIR/SAM_strand_separator.pl $SAM_file $SS_lib_type";
			&process_cmd($cmd);
		}

		push (@sam_info, [$plus_strand_sam, '+'], [$minus_strand_sam, '-']);
	}
	else {
		push (@sam_info, [$SAM_file, '+']);
	}
			

    my $thread_helper = new Thread_helper($CPU);
    
	foreach my $sam_info_aref (@sam_info) {
				
		my ($sam, $strand) = @$sam_info_aref;
				
		my $thread = threads->create('prep_read_partitions', $sam, $strand);
        #push (@threads, $thread);

        $thread_helper->add_thread($thread);
        
	}
    
    $thread_helper->wait_for_all_threads_to_complete();

    my @failures = $thread_helper->get_failed_threads();
    my $ret = 0;
    if (@failures) {
        foreach my $thread (@failures) {
            if (my $error = $thread->error()) {
                print STDERR "Error, thread exited with error $error\n";
                $ret++;
            }
        }
    }
    
	print "##\nDone\n##\n\n" unless($ret);
    
	exit($ret);

	
	

}


sub prep_read_partitions {
    my ($sam, $strand) = @_;

    ## define fragments
    my $cmd = "$UTIL_DIR/SAM_to_frag_coords.pl --CPU $CPU --sort_buffer $sort_buffer --sam $sam --min_insert_size 1 --max_insert_size $max_intron_length "; ## writes file: $sam_file.frag_coords
    &process_cmd($cmd) unless (-s "$sam.frag_coords");
    
    ## define coverage
    $cmd = "$UTIL_DIR/fragment_coverage_writer.pl $sam.frag_coords > $sam.frag_coverage.wig";

    unless (-s "$sam.frag_coverage.wig.ok") {
        &process_cmd($cmd);
        &process_cmd("touch $sam.frag_coverage.wig.ok");
    }

    my $partitions_file = "$sam.minC$min_coverage.gff";

    ## define partitions
    $cmd = "$UTIL_DIR/define_coverage_partitions.pl $sam.frag_coverage.wig $min_coverage $strand > $partitions_file";
    unless (-s "$partitions_file.ok") {
        &process_cmd($cmd);
        &process_cmd("touch $partitions_file.ok");
    }


    ## extract reads per partition
    $cmd = "$UTIL_DIR/extract_reads_per_partition.pl --partitions_gff $partitions_file "
        . " --coord_sorted_SAM $sam"
        . " --parts_per_directory $parts_per_dir"
        . " --min_reads_per_partition $min_reads_per_partition ";

    if ($SS_lib_type) {
        $cmd .= " --SS_lib_type $SS_lib_type ";
    }
    
    my $partitions_dir = "Dir_" . basename($partitions_file);
    unless (-d $partitions_dir && -e "$partitions_dir.ok") {
        &process_cmd($cmd);
        &process_cmd("touch $partitions_dir.ok");
    }
    
    return;
   
    
}


####
sub process_cmd {
	my ($cmd) = @_;

	print STDERR "CMD: $cmd\n";
	
	my $ret = system($cmd);

	if ($ret) {
	    confess "Error, command $cmd died with ret $ret";
	}

	return;
}
