#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Basename;
use Cwd;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

######################################################################
#
#  Required:
#  --genome <string>           target genome to align to
#  --transcripts <string>      cdna sequences to align
#
#  Optional:
#  -N <int>                    number of top hits (default: 1)
#  -I <int>                    max intron length
#  --CPU <int>                 number of threads (default: 2)
#  --no_chimera                do not report chimeric alignmetnts
#  --SAM                       output in SAM format
#  --gtf <string>              gene structure annotations in gtf format (for genome building)
#  --splice_assist <string>    splice assist mode (introns or splicesites)
#  
#
#######################################################################


__EOUSAGE__

    ;


my ($genome, $transcriptDB, $max_intron);
my $CPU = 2;

my $help_flag;

my $number_top_hits = 1;
my $no_chimera_flag = 0;
my $SAM_flag = 0;
my $gtf_file = 0;
my $splice_assist = "";

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'transcripts=s' => \$transcriptDB,
             'I=i' => \$max_intron,
             'CPU=i' => \$CPU,
             'N=i' => \$number_top_hits,
             'no_chimera' => \$no_chimera_flag,
             'SAM' => \$SAM_flag,
             'gtf=s' => \$gtf_file,
             'splice_assist=s' => \$splice_assist,
             );


if ($help_flag) {
    die $usage;
}


unless ($genome && $transcriptDB) {
    die $usage;
}

if ($splice_assist) {
    unless ($splice_assist =~ /^(introns|splicesites)$/) {
        die "Error, splice_assist $splice_assist not recognized as an option";
    }
}


main: {
	
	my $genomeName = basename($genome);
	my $genomeDir = $genomeName . ".gmap";

	my $genomeBaseDir = dirname($genome);

	my $cwd = cwd();
	
	unless (-d "$genomeBaseDir/$genomeDir") {
		
        #my $cmd = "gmap_build -D $genomeBaseDir -d $genomeBaseDir/$genomeDir -k 13 $genome >&2";
        my $cmd = "gmap_build -D $genomeBaseDir -d $genomeDir -k 13 $genome >&2";
		&process_cmd($cmd);
        

        if ($gtf_file) {
            &build_intron_n_splice_info_files($genomeBaseDir, $genomeDir, $gtf_file);
        }
        
	}

	
	## run GMAP

    my $num_gmap_top_hits = $number_top_hits;
    if ((! $no_chimera_flag) && $num_gmap_top_hits == 1) {
        $num_gmap_top_hits = 0; # reports two hits if chimera with this setting.
    }
    
    my $format = ($SAM_flag) ? "samse" : "3";

	my $cmd = "gmap -D $genomeBaseDir -d $genomeDir $transcriptDB -f $format -n $num_gmap_top_hits -x 50 -t $CPU -B 5 ";
	if ($max_intron) {
		$cmd .= " --intronlength=$max_intron ";
	}

    if ($splice_assist) {
        $cmd .= " -m ref_${splice_assist}.iit";
    }
	
	&process_cmd($cmd);

	exit(0);
}

####
sub process_cmd {
	my ($cmd) = @_;
	
	print STDERR "CMD: $cmd\n";
	#return;

	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret)";
	}

	return;
}

####
sub build_intron_n_splice_info_files {
    my ($genomeBaseDir, $genomeDir, $gtf_file) = @_;

    my $genome_maps_dir = "$genomeBaseDir/$genomeDir/$genomeDir.maps";

    my $introns_iit_file = "$genome_maps_dir/ref_introns.iit";

    my $cmd = "gtf_introns < $gtf_file | iit_store -o $introns_iit_file";
    &process_cmd($cmd);

    if (! -e $introns_iit_file) {
        die "Error, no introns $introns_iit_file";
    }

    my $splicesites_iit_file = "$genome_maps_dir/ref_splicesites.iit";

    $cmd = "gtf_splicesites < $gtf_file | iit_store -o $splicesites_iit_file";
    &process_cmd($cmd);

    if (! -e $splicesites_iit_file) {
        die "Error, no splicesites $splicesites_iit_file";
    }

    return;
}

