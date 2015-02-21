#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use File::Basename;
use Cwd;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $HISAT_HOME = $ENV{HISAT_HOME} or die "Error, need env var HISAT_HOME set to the HISAT installation directory.\n\n";
 


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
#  --out_prefix <string>       output prefix (default: hisat)
#
#######################################################################


__EOUSAGE__

    ;


my ($genome, $reads);

my $CPU = 2;

my $help_flag;

my $out_prefix = "hisat";
my $gtf_file;
my $no_sarray = "";

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


main: {
	
    my $hisat_index = "$genome.hisat.idx";
    if (! -s "$hisat_index.1.bt2") {
        ## build hisat index

        my $cmd = "$HISAT_HOME/hisat-build $genome $hisat_index";
        &process_cmd($cmd);
    }

    my $gtf_splice = "$gtf_file.hisat.splice";
    my $splice_incl = "";
    
    if ($gtf_file) {

        unless (-s $gtf_splice) {
            my $cmd = "$HISAT_HOME/extract_splice_sites.py $gtf_file > $gtf_file.hisat.splice";
            &process_cmd($cmd);
        }
        
        $splice_incl = " --known-splicesite-infile $gtf_splice ";
    }
    
    ## run HISAT
    
    $reads = &add_zcat_fifo_and_add_hisat_params($reads);

    my $cmd = "$HISAT_HOME/hisat -x $hisat_index -q $reads $splice_incl -p $CPU  | samtools view -@ $CPU -F 4 -Sb - | samtools sort -@ $CPU -no - - > $out_prefix.cSorted.bam";
    
    &process_cmd($cmd);

    if (-s "$out_prefix.cSorted.bam") {
        $cmd = "samtools index $out_prefix.cSorted.bam";
        &process_cmd($cmd);
    }
    
	exit(0);
}


####
sub add_zcat_fifo_and_add_hisat_params {
    my ($reads) = @_;

    $reads =~ s/^\s+|\s+$//g;
    
    my @adj_reads_list;

    my $counter = 0;
    my @read_files = split(/\s+/, $reads);
    
    foreach my $reads_file (@read_files) {
        
        $counter++;
        
        if ($reads_file =~ /\.gz$/) {
            $reads_file = "<(zcat $reads_file)";
        }
       
        # add decoration 
        $reads_file = (scalar(@read_files) == 2) ? "-$counter $reads_file" : "-U $reads_file";

        push (@adj_reads_list, $reads_file);
    }
    
    my $adj_reads = join(" ", @adj_reads_list);

    return($adj_reads);
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



