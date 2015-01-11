#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Cwd;
use File::Basename;
use Carp;
use Data::Dumper;

use Getopt::Long qw(:config no_ignore_case bundling);

$ENV{LC_ALL} = 'C';  # critical for proper sorting using [system "sort -k1,1 ..."] within the perl script

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  --left and --right <string>   (if paired reads)
#     or
#  --single <string>             (if unpaired reads)
#
#  Required inputs:
#   
#  --target <string>            multi-fasta file containing the target sequences (should be named {refName}.fa )
#
#  --out_prefix|o               output prefix (default: bwa)
#
# 
#  ## General options
#
#  Any options after '--' are passed onward to the alignments programs (except BLAT -which has certain options exposed above).
#
#  To set the number of processors used by BWA use:
#     -- -t 16
#  You could also set other options such as 'mismatch penalty' (via -- -M INT), etc. 
#
####################################################################################################################



_EOUSAGE_

    ;


my $help_flag;
my $target_db; 
my $left_file;
my $right_file;
my $single_file;

my $output_prefix = "bwa";


unless (@ARGV) {
    die $usage;
}

&GetOptions ( 'h' => \$help_flag,
              
              ## required inputs
              'left=s' => \$left_file,
              'right=s' => \$right_file,
              
              'single=s' => \$single_file,
              

              "target=s" => \$target_db,

              'output_prefix|o=s' => \$output_prefix,
              
    );






if ($help_flag) { die $usage; }

unless ($target_db && -s $target_db) { 
    die $usage . "Must specify target_db and it must exist at that location";
}


unless ( ($single_file && -e $single_file)
         || 
         ($left_file && -e $left_file
          && $right_file && -e $right_file)) {
    die $usage . "sorry, cannot find $left_file and $right_file";
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";

    my $ret = system($cmd);

    if ($ret) {
        confess "Error, cmd: $cmd died with ret $ret";
    }

    return;
}




main: {

    
    my $cmd = "bwa index $target_db";
    &process_cmd($cmd) unless (-e "$target_db.bwt");;
    
    $cmd = "samtools faidx $target_db";
    &process_cmd($cmd) unless (-e "$target_db.fai");
    

    if ($left_file && $right_file) {

        $cmd = "bwa aln @ARGV $target_db $left_file > $left_file.sai";
        &process_cmd($cmd);
        
        $cmd = "bwa aln @ARGV $target_db $right_file > $right_file.sai";
        &process_cmd($cmd);
        
        $cmd = "bwa sampe $target_db $left_file.sai $right_file.sai $left_file $right_file | samtools view -bS -F 4 - | samtools sort - $output_prefix";
        &process_cmd($cmd);
        
    }
    else {
                
        $cmd = "bwa aln @ARGV $target_db $single_file > $single_file.sai";
        &process_cmd($cmd);
        
        $cmd = "bwa samse $target_db $single_file.sai $single_file | samtools view -bS -F 4 - | samtools sort - $output_prefix";
        &process_cmd($cmd);
    }
    
    
    exit(0);
}
    
