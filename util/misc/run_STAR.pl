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
#
#######################################################################


__EOUSAGE__

    ;


my ($genome, $reads);

my $CPU = 2;

my $help_flag;

my $out_prefix = "star";
my $gtf_file;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'reads=s' => \$reads,
             'CPU=i' => \$CPU,
             'out_prefix=s' => \$out_prefix,
             'G=s' => \$gtf_file,

    );


unless ($genome && $reads) {
    die $usage;
}


my $star_prog = `which STAR`;
chomp $star_prog;
unless ($star_prog =~ /\w/) {
    die "Error, cannot locate STAR program. Be sure it's in your PATH setting.  ";
}


main: {
	
    my $star_index = "$genome.star.idx";
    if (! -s "$star_index/genomeParameters.txt") {
        ## build star index
        unless (-d $star_index) {
            mkdir($star_index) or die "Error, cannot mkdir $star_index";
        }
        
        
        my $cmd = "$star_prog --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index "
            . " --genomeFastaFiles $genome --sjdbOverhang 100 ";
        if ($gtf_file) {
            $cmd .= " --sjdbGTFfile $gtf_file ";
        }
        
        &process_cmd($cmd);
    }

        
    ## run STAR
    
    my @tmpfiles;
    
    my $pipeliner = new Pipeliner(-verbose => 1);

    my $cmd = "$star_prog --runThreadN $CPU --genomeDir $star_index --outSAMtype BAM SortedByCoordinate --readFilesIn $reads 
";
    if ($reads =~ /\.gz$/) {
        $cmd .= " --readFilesCommand 'gzip -c' ";
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



