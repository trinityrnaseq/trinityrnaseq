#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use File::Basename;
use FindBin;
use Cwd;

my $usage = <<__EOUSAGE__;

###############################################################################
#
# Required:
#
#  --target <string>       : genome multifasta file
#  --seqType <string>      :type of reads ('fq' or 'fa')
#
#  If paired reads:
#      --left  <string>    :left reads
#      --right <string>    :right reads
#
#  Or, if unpaired reads:
#      --single <string>   :single reads 
#
#  --output|-o <string>               :name of directory for output 
#
#  -i <int>                :min intron length
#  -I <int>                :max intron length 
#
#
# Optional:
#
#  --SS_lib_type <string>          :strand-specific library type : {RF, FR, R, F}  (RF =~ fr-firststrand)
#  --CPU <int>                     :number of CPUs to use
#
#  --GTF <string>                  :annotations to assist, provided in GTF file format
#
###############################################################################


__EOUSAGE__

    ;



my $help_flag;

my ($seqType, $left_file, $right_file, $single_file, $SS_lib_type, $output_dir, $CPU, $target_genome, $min_intron_length, $max_intron_length);

my $GTF_annots;

&GetOptions ( 'h' => \$help_flag,
              
              ## general opts
              "seqType=s" => \$seqType,
              "left=s" => \$left_file,
              "right=s" => \$right_file,
              "single=s" => \$single_file,
              "target=s" => \$target_genome,

              "SS_lib_type=s" => \$SS_lib_type,
              "output|o=s" => \$output_dir,

              'CPU=i' => \$CPU,
              'GTF=s' => \$GTF_annots,

              'i=i' => \$min_intron_length,
              'I=i' => \$max_intron_length,
              
              );


if (@ARGV) {
    die "Error, do not understand options: @ARGV\n";
}

unless ($seqType && $target_genome && ( ($single_file) || ($left_file && $right_file) ) && $output_dir) {
    die $usage;
}

if ($GTF_annots) {
    ## add full path
    if ($GTF_annots !~ /^\//) {
        $GTF_annots = cwd() . "/$GTF_annots";
    }
}


main: {

    my $util_dir = "$FindBin::RealBin/..";
    
    #############################
    ## align reads using Tophat
    #############################
    
    my $cmd = "$util_dir/alignReads.pl --seqType $seqType --output $output_dir --aligner tophat2 --target $target_genome ";
    if ($single_file) {
        $cmd .= " --single $single_file ";
    }
    else {
        $cmd .= " --left $left_file --right $right_file ";
    }
    if ($SS_lib_type) {
        $cmd .= " --SS_lib_type $SS_lib_type ";
    }
    if ($min_intron_length) {
        $cmd .= " -i $min_intron_length ";
    }
    if ($max_intron_length) {
        $cmd .= " -I $max_intron_length ";
    }
    

    if ($CPU || $GTF_annots) {
        $cmd .= " -- ";

        if ($CPU) {
            $cmd .= " -p $CPU ";
        }
        if ($GTF_annots) {
            $cmd .= " --GTF $GTF_annots ";
            
            my $transcriptome_index_dir = "$GTF_annots.tuxedo_indices";
            if (-d $transcriptome_index_dir) {
                my $prefix = basename($GTF_annots);
                $prefix =~ s/\.[^\.]+$//;
                $cmd .= " --transcriptome-index $transcriptome_index_dir/$prefix ";
            }
            else {
                #die "Error, cannot find: $transcriptome_index_dir";
                print STDERR "WARNING: cannot find pre-built index for $GTF_annots; this will slow it down.\n";
            }
                        
        }
    }
    
    &process_cmd($cmd);
    
    
    #######################################
    # assemble transcripts using cufflinks
    #######################################

    my $tophat_alignment_bam = "$output_dir/tophat_out/accepted_hits.bam";

    $cmd = "cufflinks  --no-update-check --upper-quartile-norm -o $output_dir/cuff_out ";
    if ($GTF_annots) {
        $cmd .= " --GTF-guide $GTF_annots ";
    }


    
    
    if ($SS_lib_type) {
        my %tuxedo_lib_type = (RF => "fr-firststrand",
                               R => "fr-firststrand",
                               FR => "fr-secondstrand",
                               F => "fr-secondstrand",
                               );
        
        my $lib_type = $tuxedo_lib_type{$SS_lib_type} or die "Error, cannot determine tophat lib_type based on SS_lib_type: $SS_lib_type ";
        
        $cmd .= " --library-type $lib_type ";
    }
    
    if ($CPU) {
        $cmd .= " -p $CPU ";
    }
    
    $cmd .= " $tophat_alignment_bam ";
    
    &process_cmd($cmd);

    exit(0);
    
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}



