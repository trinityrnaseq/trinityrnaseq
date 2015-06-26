#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::Bin/../../PerlLib");
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
#  --reads  <string>           fastq files. If pairs, indicate both in quotes, ie. "left.fq right.fq"
#
#  Optional:
#  -G <string>                 GTF file for incorporating reference splice site info.
#  --CPU <int>                 number of threads (default: 2)
#  --out_prefix <string>       output prefix (default: star)
#  --out_dir <string>          output directory (default: current working directory)
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

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'reads=s' => \$reads,
             'CPU=i' => \$CPU,
             'out_prefix=s' => \$out_prefix,
             'G=s' => \$gtf_file,
             'out_dir=s' => \$out_dir,
    );


unless ($genome && $reads) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}


my $star_prog = `which STAR`;
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
            . " --twopassMode Basic "
            . " --genomeFastaFiles $genome --sjdbOverhang 100 ";
        if ($gtf_file) {
            $cmd .= " --sjdbGTFfile $gtf_file ";
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
        . " --readFilesIn $reads "
        . " --chimJunctionOverhangMin 12 "
        . " --chimSegmentMin 12 "
        . " --outSAMstrandField intronMotif "
        . " --twopassMode Basic "
        . " --alignSJDBoverhangMin 10 ";
    
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



