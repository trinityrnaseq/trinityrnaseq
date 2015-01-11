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
#  -I <int>                    max intron length
#  --CPU <int>                 number of threads (default: 2)
#
#######################################################################


__EOUSAGE__

    ;


my ($genome, $transcriptDB, $max_intron);
my $CPU = 2;

my $help_flag;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'transcripts=s' => \$transcriptDB,
             'I=i' => \$max_intron,
             'CPU=i' => \$CPU,
             );


unless ($genome && $transcriptDB) {
    die $usage;
}


main: {
	
	my $genomeName = basename($genome);
	my $genomeDir = $genomeName . ".gmap";

	my $genomeBaseDir = dirname($genome);

	my $cwd = cwd();
	
	unless (-d "$genomeBaseDir/$genomeDir") {
		
        my $cmd = "gmap_build -D $genomeBaseDir -d $genomeDir -k 13 $genome >&2";
		&process_cmd($cmd);
	}

	
	## run GMAP
	my $cmd = "gmap -D $genomeBaseDir -d $genomeDir $transcriptDB -f 3 -n 0 -x 50 -t $CPU ";
	if ($max_intron) {
		$cmd .= " --intronlength=$max_intron ";
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



