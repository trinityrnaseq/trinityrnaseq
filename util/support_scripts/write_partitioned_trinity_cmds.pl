#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

####################################################################################
#
#  usage: $0 --reads_list_file <string> [Trinity params]
#
# Required:
#
# --reads_list_file <string>      file containing list of filenames corresponding 
#                                  to the reads.fasta
#
# Optional for singularity usage:
#
#   --singularity_img <string>     path to singularity image
#
#   --singularity_extra_params <string>  singularity extra parameters to include
#
#####################################################################################


__EOUSAGE__

    ;


my $reads_file;
my $help_flag;
my $singularity_img;
my $singularity_extra_params = "";


&GetOptions (
    'reads_list_file=s' => \$reads_file,
    'h' => \$help_flag,
    'singularity_img=s' => \$singularity_img,
    'singularity_extra_params=s' => \$singularity_extra_params,    
    );


my @TRIN_ARGS = @ARGV;

if ($help_flag) {
    die $usage;
}

unless ($reads_file && -e $reads_file) {
    die $usage;
}

unless (-s $reads_file) {
    die "Error, reads file listing: $reads_file is empty.  This tends to happen when there were too few reads to assemble. ";
}


my $trin_args = "";
while (@TRIN_ARGS) {
    my $arg = shift @TRIN_ARGS;
    
    if ($arg =~ /bfly_opts/) {
        my $val = shift @TRIN_ARGS;
        # retain quotes around multiparams
        
        $trin_args .= "$arg \"$val\" ";
    }
    else {
        $trin_args .= "$arg ";
    }
}


open (my $fh, $reads_file) or die "Error, cannot open file $reads_file";
while (<$fh>) {
	chomp;
    my @x = split(/\s+/);
    
    my $file = pop @x;
    

    my $cmd = "$FindBin::RealBin/../../Trinity --single \"$file\" --output \"$file.out\" $trin_args ";

    if ($singularity_img) {
        $cmd = "singularity exec $singularity_extra_params $singularity_img $cmd";
    }
    print "$cmd\n";
}

exit(0);




		
