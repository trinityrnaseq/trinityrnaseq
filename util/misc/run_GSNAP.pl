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
#  --reads  <string>           fastq files. If pairs, indicate both in quotes, ie. "left.fq right.fq"
#
#  Optional:
#  -N <int>                    number of top hits (default: 1)
#  -I <int>                    max intron length (default: 500000)
#  --CPU <int>                 number of threads (default: 2)
#  --out_prefix <string>       output prefix (default: gsnap)
#
#######################################################################


__EOUSAGE__

    ;


my ($genome, $reads);

my $max_intron = 500000;
my $CPU = 2;

my $help_flag;

my $num_top_hits = 1;
my $out_prefix = "gsnap";

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'reads=s' => \$reads,
             'I=i' => \$max_intron,
             'CPU=i' => \$CPU,
             'N=i' => \$num_top_hits,
             'out_prefix=s' => \$out_prefix,
             );


unless ($genome && $reads) {
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
    
    my $cmd = "gsnap -D $genomeBaseDir -d $genomeDir -A sam -N 1 -w $max_intron -n $num_top_hits -t $CPU $reads | samtools view -bS - | samtools sort -@ $CPU - $out_prefix.cSorted";
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



