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
#  --gtf <string>              gene structure annotations in gtf format
#  --CPU <int>                 number of threads (default: 2)
#  -o|--output <string>        bam output filename (default: basename(transcripts).mm2.bam)
#
#######################################################################


__EOUSAGE__

    ;

my $genome;
my $transcripts;
my $gtf;
my $CPU = 2;

my $help_flag;
my $output;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'transcripts=s' => \$transcripts,
             'gtf=s' => \$gtf,
             'CPU=i' => \$CPU,
             'o|output=s' => \$output,
    );


if ($help_flag) {
    die $usage;
}


unless ($genome && $transcripts) {
    die $usage;
}

unless ($output) {
    $output = basename($transcripts) . ".mm2.bam";
}

main: {

	
	my $genomeBaseDir = dirname($genome);
	my $genomeName = basename($genome);
	my $genomeDir = "$genomeBaseDir/$genomeName" . ".mm2";
    

    unless (-d $genomeDir) {
        &process_cmd("mkdir $genomeDir");
    }
    
	my $cwd = cwd();
	
    my $mm2_idx = "$genomeDir/$genomeName.mmi";
	
    my $splice_file = "";
    if ($gtf) {
        $splice_file = "$genomeDir/anno.bed";
    }


    unless (-e $mm2_idx) {
        
        my $cmd = "minimap2 -d $mm2_idx $genome";
        &process_cmd($cmd);

        if ($gtf) {
            my $cmd = "paftools.js gff2bed $gtf > $splice_file";
            &process_cmd($cmd);
        }
	
    }
    
	
	## run minimap2

    my $splice_param = "";
    if ($splice_file) {
        $splice_param = "--junc-bed $splice_file";
    }

    my $cmd = "minimap2 -ax splice $splice_param --secondary=no -O6,24 -B4 -L -t $CPU -cs -ub $mm2_idx $transcripts | samtools view -Sb | samtools sort -o $output && samtools index $output";
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



