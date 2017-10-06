#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::RealBin/../../PerlLib");
use Pipeliner;
use File::Basename;
use Cwd;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

######################################################################
#
#  Required:
#  --genome <string>           target genome to align to
#  --gtf <string>              annotations in gtf format
#  --samples_file  <string>    trinity samples file
#  
#  Optional:
#  --CPU <int>                 number of threads (default: 2)
#
#######################################################################


__EOUSAGE__

    ;


my ($genome);
my $samples_file;

my $CPU = 2;

my $help_flag;
my $gtf_file;


&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'samples_file=s' => \$samples_file,
             'CPU=i' => \$CPU,
             'gtf=s' => \$gtf_file,
    );


unless ($genome && $samples_file && $gtf_file) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, cannot recognize opts: @ARGV";
}


my $star_prog = `which STAR`;
chomp $star_prog;
unless ($star_prog =~ /\w/) {
    die "Error, cannot locate STAR program. Be sure it's in your PATH setting.  ";
}


main: {
	
    ## ensure all full paths
    $genome = &Pipeliner::ensure_full_path($genome);
    $gtf_file = &Pipeliner::ensure_full_path($gtf_file);

        
    my @read_sets = &parse_samples_file($samples_file);    

    my $pipeliner = new Pipeliner(-verbose => 1);
    my $star_index = "$genome.star.idx";
    my $star_index_chkpt = "$star_index/build.ok";
    ## build star index
    unless (-d $star_index) {
        mkdir($star_index) or die "Error, cannot mkdir $star_index";
    }
    
    my $cmd = "$star_prog --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index "
        . " --genomeFastaFiles $genome "
        . " --limitGenomeGenerateRAM 40419136213 "
        . " --sjdbGTFfile $gtf_file "
        . " --sjdbOverhang 150 ";

    $pipeliner->add_commands( new Command($cmd, $star_index_chkpt));

    $pipeliner->run();
    
    my $checkpoint_dir = "star_aln_chkpts." . basename($genome);
    unless (-d $checkpoint_dir) {
        mkdir($checkpoint_dir) or die "Error, cannot mkdir $checkpoint_dir";
    }
            
    foreach my $read_set_aref (@read_sets) {
        my ($sample_id, $left_fq, $right_fq) = @$read_set_aref;
            
        my $cmd = "$star_prog "
            . " --runThreadN $CPU "
            . " --genomeDir $star_index "
            . " --outSAMtype BAM SortedByCoordinate "
            . " --runMode alignReads "
            . " --readFilesIn $left_fq $right_fq "
            . " --twopassMode Basic "
            . " --alignSJDBoverhangMin 10 "
            . " --limitBAMsortRAM 20000000000";
        
            
        if ($left_fq =~ /\.gz$/) {
            $cmd .= " --readFilesCommand 'gunzip -c' ";
        }
        
        $pipeliner->add_commands( new Command($cmd, "$checkpoint_dir/star_align.$sample_id." . basename($genome) . ".ok") );
        
        my $bam_outfile = "Aligned.sortedByCoord.out.bam";
        my $renamed_bam_outfile = "$sample_id.cSorted.star." . basename($genome) . ".bam";
        $pipeliner->add_commands( new Command("mv $bam_outfile $renamed_bam_outfile", "$checkpoint_dir/$renamed_bam_outfile.ok") );
        
        $pipeliner->add_commands( new Command("samtools index $renamed_bam_outfile", "$checkpoint_dir/$renamed_bam_outfile.bai.ok") );
    
    
        $pipeliner->run();
    }
    
    
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
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples;

    open(my $fh, $samples_file) or die "Error, cannot open file $samples_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($cond, $rep, $fq_a, $fq_b) = @x;

        push (@samples, [$rep, $fq_a, $fq_b]);
    }
    close $fh;

    return (@samples);
}

