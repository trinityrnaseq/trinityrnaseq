#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::RealBin/../../PerlLib");
use Pipeliner;
use File::Basename;
use Cwd;
use List::Util qw(min);

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

######################################################################
#
#  Required:
#  --bam_file_list  <string>   bam_file_list (format: sample_id(tab)/path/to/file.bam)
#  
#  Optional
#  --gtf <string>              annotations in gtf format
#  --CPU <int>                 number of threads (default: 2)
#  --SS_lib_type <string>      
# 
#
#######################################################################


__EOUSAGE__

    ;


my $bam_file_list;

my $CPU = 2;

my $help_flag;
my $gtf_file;

my $SS_lib_type;

&GetOptions( 'h' => \$help_flag,

             'bam_file_list=s' => \$bam_file_list,
             'CPU=i' => \$CPU,
             'gtf=s' => \$gtf_file,

             'SS_lib_type=s' => \$SS_lib_type,
    );


unless ($bam_file_list) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, cannot recognize opts: @ARGV";
}



main: {
	
    ## ensure all full paths
    $gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;
    
    my @entries;
    open(my $fh, $bam_file_list) or die $!;
    while(<$fh>) {
        chomp;
        my ($sample_id, $bam_filename) = split(/\t/);
        unless ($sample_id && $bam_filename) {
            die "Error, bam file list doesn't have expected tab-delimited formatting at $_";
        }
        unless (-e $bam_filename) {
            die "Error, cannot locate file: $bam_filename";
        }
        
        push (@entries, { sample_id => $sample_id,
                          bam => $bam_filename } );

    }

    my $pipeliner = new Pipeliner(-verbose => 1);
      
    
    my @indiv_stringtie_outputs;
    
    foreach my $entry (@entries) {
        my $cmd = "stringtie " . $entry->{bam} . 
            " -l STRG." . $entry->{sample_id} . " ";
        
        if ($SS_lib_type) {
            if ($SS_lib_type =~ /^R/i) {
                $cmd .= " --rf ";
            }
            elsif ($SS_lib_type =~ /^F/i) {
                $cmd .= " --fr ";
            }
            else {
                die "Error, cannot determine strand-specificity type from $SS_lib_type";
            }
        }

        $cmd .= " -o " . $entry->{sample_id} . ".strg.gtf";
        
        push (@indiv_stringtie_outputs, $entry->{sample_id} . ".strg.gtf");

        if ($gtf_file) {
            $cmd .= " -G $gtf_file ";
        }
        my $checkpoint = $entry->{sample_id} . ".strg.ok";
        
        $pipeliner->add_commands( new Command($cmd, $checkpoint));

    }
    
    $pipeliner->run();


    # merge transcripts
    if (scalar(@indiv_stringtie_outputs) > 1) {
        my $cmd = "stringtie --merge -o stringtie.merged.gtf -g 1 ";
        if ($gtf_file) {
            $cmd .= " -G $gtf_file ";
        }
        $cmd .= join(" ", @indiv_stringtie_outputs);
    
        $pipeliner->add_commands( new Command($cmd, "stringtie_merge.ok"));
        
        $pipeliner->run();
        
    }
    
        
	exit(0);
}


