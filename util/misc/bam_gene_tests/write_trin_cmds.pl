#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = <<__EOUSAGE__;

####################################################################################
#
#  usage: $0 --reads_list_file <string> --out_token <string> [Trinity params]
#
# Required:
#
# --reads_list_file <string>      file containing list of filenames corresponding 
#                                  to the reads.fasta
#
# --out_token <string>            token added to the output file name
#
#####################################################################################

#  Example: 
#
#    write_trin_cmds.pl  --reads_list_file ReadPartitions.listing --out_token origbfly --SS_lib_type F --full_cleanup_ET --CPU 1 --bfly_jar ~/SVN/trinityrnaseq/trunk/Butterfly/Butterfly.jar --JM 1G --seqType fa


__EOUSAGE__

    ;


my $reads_file;
my $help_flag;
my $out_token = "";


&GetOptions (
             
             'reads_list_file=s' => \$reads_file,
             'h' => \$help_flag,
             
             'out_token=s' => \$out_token,
             
             
             );

my @TRIN_ARGS = @ARGV;

if ($help_flag) {
    die $usage;
}

unless ($reads_file && -s $reads_file) {
    die $usage;
}

unless ($out_token) {
    die $usage;
}

my $trin_args = join(" ", @TRIN_ARGS);


open (my $fh, $reads_file) or die "Error, cannot open file $reads_file";
while (<$fh>) {
	chomp;
    my @x = split(/\s+/);
    
    my $file = pop @x;
    
    my $cmd = "$FindBin::Bin/../../../Trinity --single \"$file\" --output \"$file.trinity.$out_token\" $trin_args ";
    
    print "$cmd\n";
}

exit(0);




		
