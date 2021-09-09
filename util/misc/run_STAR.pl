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
#     --gtf|G <string>                 GTF file for incorporating reference splice site info. (recommended)
#
#  --reads  <string>           fastq files. If pairs, indicate both in quotes, ie. "left.fq right.fq"
#
#  Optional:
#  --CPU <int>                 number of threads (default: 2)
#  --out_prefix <string>       output prefix (default: star)
#  --out_dir <string>          output directory (default: current working directory)
#  --star_path <string>        full path to the STAR program to use.
#  --patch <string>            genomic targets to patch the genome fasta with.
#  --chim_search               include Chimeric.junction outputs
#  --max_intron <int>          max intron length (and PE gap size)
#
#######################################################################


__EOUSAGE__

    ;


my ($genome, $reads);

my $CPU = 2;

my $help_flag;

my $out_prefix = "star";
my $gtf_file;
my $out_dir;
my $ADV = 0;

my $star_path = "STAR";
my $patch;
my $chim_search;
my $max_intron;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'reads=s' => \$reads,
             'CPU=i' => \$CPU,
             'out_prefix=s' => \$out_prefix,
             'gtf|G=s' => \$gtf_file,
             'out_dir=s' => \$out_dir,
             'ADV' => \$ADV,
             'star_path=s' => \$star_path,
             'patch=s' => \$patch,
             'chim_search' => \$chim_search,
             "max_intron=i" => \$max_intron,
    );


unless ($genome && $reads) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, cannot recognize opts: @ARGV";
}


my $star_prog = `which $star_path`;
chomp $star_prog;
unless ($star_prog =~ /\w/) {
    die "Error, cannot locate STAR program. Be sure it's in your PATH setting.  ";
}


main: {
	
    ## ensure all full paths
    $genome = &Pipeliner::ensure_full_path($genome);
    $gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;

    my @read_files = split(/\s+/, $reads);
    foreach my $read_file (@read_files) {
        if ($read_file) {
            $read_file = &Pipeliner::ensure_full_path($read_file);
        }
    }
    $reads = join(" ", @read_files);
    
    if ($out_dir) {
        unless (-d $out_dir) {
            mkdir $out_dir or die "Error, cannot mkdir $out_dir";
        }
        chdir $out_dir or die "Error, cannot cd to $out_dir";
    }
    

    my $star_index = "$genome.star.idx";
    if (! -e "$star_index/build.ok") {
        ## build star index
        unless (-d $star_index) {
            mkdir($star_index) or die "Error, cannot mkdir $star_index";
        }
        
        
        my $cmd = "$star_prog --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index "
            . " --genomeFastaFiles $genome "
            . " --limitGenomeGenerateRAM 40419136213 ";
        if ($gtf_file) {
            $cmd .= " --sjdbGTFfile $gtf_file "
                . " --sjdbOverhang 100 ";
            
        }
        
        &process_cmd($cmd);
        
        &process_cmd("touch $star_index/build.ok");
        
    }

        
    ## run STAR
    
    my @tmpfiles;
    
    my $pipeliner = new Pipeliner(-verbose => 1);

    my $cmd = "$star_prog "
        . " --runThreadN $CPU "
        . " --genomeDir $star_index "
        . " --outSAMtype BAM SortedByCoordinate "
        . " --runMode alignReads "
        . " --readFilesIn $reads "
        . " --twopassMode Basic "
        . " --alignSJDBoverhangMin 10 "
        . " --outSAMstrandField intronMotif "
        . " --outSAMunmapped Within "
        . " --alignInsertionFlush Right "
        . " --alignSplicedMateMapLminOverLmate 0 "
        . " --alignSplicedMateMapLmin 30 "
        . " --alignSJstitchMismatchNmax 5 -1 5 5 "  #which allows for up to 5 mismatches for non-canonical GC/AG, and AT/AC junctions, and any number of mismatches for canonical junctions (the default values 0 -1 0 0 replicate the old behavior (from AlexD)      
        . " --peOverlapNbasesMin 12 "
        . " --peOverlapMMp 0.1 "
        . " --limitBAMsortRAM 20000000000";


    if ($max_intron) {
    
        $cmd .= " --alignMatesGapMax $max_intron "
            . " --alignIntronMax $max_intron ";
    }
    
    
    if ($chim_search) {
        $cmd .= " --chimJunctionOverhangMin 8 "
             .  " --chimOutJunctionFormat 1 "
             .  " --chimSegmentMin 12 "
             .  " --chimSegmentReadGapMax parameter 3 "
             .  " --chimMultimapNmax 20 "
             .  " --chimOutType Junctions WithinBAM "
             .  " --chimScoreJunctionNonGTAG -4 "
             .  " --chimNonchimScoreDropMin 10 "
             .  " --chimMultimapScoreRange 3 ";
    }
    
    if ($patch) {
        $cmd .= " --genomeFastaFiles $patch ";
    }
        
        
    
    if ($reads =~ /\.gz$/) {
        $cmd .= " --readFilesCommand 'gunzip -c' ";
    }

    $pipeliner->add_commands( new Command($cmd, "star_align.ok") );
    


    my $bam_outfile = "Aligned.sortedByCoord.out.bam";
    my $renamed_bam_outfile = "$out_prefix.sortedByCoord.out.bam";
    $pipeliner->add_commands( new Command("mv $bam_outfile $renamed_bam_outfile", "$renamed_bam_outfile.ok") );
    
    
    $pipeliner->add_commands( new Command("samtools index $renamed_bam_outfile", "$renamed_bam_outfile.bai.ok") );
    
    
    $pipeliner->run();

    
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



