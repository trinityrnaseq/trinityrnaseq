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
#  -I|--max_intron_length <int>   maximum intron length (default: 100000)
#
#  --incl_out_gff3             include gff3 formatted output file for alignments.
#  --allow_secondary           allow secondary alignments (default secondary=no)
#
#  --eqx                       include --eqx flag
#  --cs                        include long format via --cs w/ minimap2
#  --hq                        pacbio CCS reads (--splice:hq for mm2)
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
my $max_intron_length = 100000;
my $incl_out_gff3;
my $allow_secondary = 0;
my $include_cs_flag = 0;
my $include_eqx_flag = 0;
my $include_hq_flag = 0;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'transcripts=s' => \$transcripts,
             'gtf=s' => \$gtf,
             'CPU=i' => \$CPU,
             'o|output=s' => \$output,
             'I|max_intron_length=i' => \$max_intron_length,
             'incl_out_gff3' => \$incl_out_gff3,
             'allow_secondary' => \$allow_secondary,
             'cs' => \$include_cs_flag,
             'hq' => \$include_hq_flag,
    );


if ($help_flag) {
    die $usage;
}


unless ($genome && $transcripts) {
    die $usage;
}


my $include_cs_param = "";
my $include_cs_token = "";
if ($include_cs_flag) {
    $include_cs_param = "--cs";
    $include_cs_token = ".cs";
}

my $include_eqx_param = "";
my $include_eqx_token = "";
if ($include_eqx_flag) {
    $include_eqx_param = "--eqx";
    $include_eqx_token = ".eqx";
}


my $include_hq_param = "";
my $include_hq_token = "";
if ($include_hq_flag) {
    $include_hq_param = ":hq";
    $include_hq_token = ".hq";
}


unless ($output) {
    $output = basename($transcripts) . ".mm2${include_cs_token}${include_eqx_token}${include_hq_token}.bam";
}




main: {

	
    my $genomeBaseDir = dirname($genome);
    my $genomeName = basename($genome);
    my $mm2_idx = "$genomeBaseDir/$genomeName" . ".mm2";
    
    my $cwd = cwd();

    my $splice_file = "$mm2_idx.splice.bed";

    unless (-e $mm2_idx) {
        
        my $cmd = "minimap2 -d $mm2_idx $genome";
        &process_cmd($cmd);
    }

    
    if ($gtf && ! -s $splice_file) {
	my $cmd = "paftools.js gff2bed $gtf > $splice_file";
	&process_cmd($cmd);
    }
	
	
    ## run minimap2

    my $splice_param = "";
    if ($splice_file) {
        $splice_param = "--junc-bed $splice_file";
    }

    my $secondary = ($allow_secondary) ? "" : "--secondary=no";

    my $cmd = "minimap2 --sam-hit-only  -ax splice$include_hq_param $splice_param $secondary -t $CPU -u b -G $max_intron_length $include_cs_param $include_eqx_param $mm2_idx $transcripts > $output.tmp.sam";
    &process_cmd($cmd);

    $cmd = "samtools view -Sb -T $genome $output.tmp.sam -o $output.tmp.unsorted.bam";
    &process_cmd($cmd);

    $cmd = "samtools sort $output.tmp.unsorted.bam -o $output";
    &process_cmd($cmd);

    $cmd = "samtools index $output";
    #&process_cmd($cmd);
    `$cmd`; # ignore error that occurs if file is too big.

    if ($incl_out_gff3) {
        $cmd = "$FindBin::Bin/SAM_to_gff3.minimap2.pl $output > $output.gff3";
        &process_cmd($cmd);
    }

    unlink("$output.tmp.sam", "$output.tmp.unsorted.bam");
    
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



