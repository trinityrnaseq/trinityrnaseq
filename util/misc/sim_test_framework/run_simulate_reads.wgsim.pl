#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$ENV{TRINITY_HOME}/PerlLib/");
use Fasta_reader;
use Cwd;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use List::Util qw (shuffle);

my $help_flag;


my $usage = <<__EOUSAGE__;

################################################################################

$0

################################################################################
#
#  * Required:
#
#  --ref_trans|R <string>                 reference transcriptome
#
#  --out_dir|O <string>                   output directory name 
#
#  --read_length <int>                    default: 76
#
#  --frag_length <int>                    default: 300
#
#  --depth_of_cov <int>                   default: 100
#
#
####
#
#  following wgsim options are pass-through:
#
# Options: 
#    -e FLOAT      base error rate [0.020]
#    -s INT        standard deviation [50]
#    -r FLOAT      rate of mutations [0.0010]
#    -R FLOAT      fraction of indels [0.15]
#    -X FLOAT      probability an indel is extended [0.30]
#    -S INT        seed for random generator [-1]
#    -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [0.05]
#    -h            haplotype mode
#    -Z INT        strand specific mode: 1=FR, 2=RF
#    -D            debug mode... highly verbose
#
#
############################################################################################


__EOUSAGE__

    ;


my $OUT_DIR;
my $ref_trans_fa;
my $read_length = 76;
my $frag_length = 300;
my $depth_of_cov = 100;


&GetOptions ( 'help' => \$help_flag,
              
              # required
              'ref_trans|R=s' => \$ref_trans_fa,
              
              # optional
              'out_dir|O=s' => \$OUT_DIR,

              'read_length=i' => \$read_length,
              'frag_length=i' => \$frag_length,
              'depth_of_cov=i' => \$depth_of_cov,
                            
    );


if ($help_flag) {
    die $usage;
}


unless ($ref_trans_fa && $OUT_DIR) {
    die $usage;
}



unless ($ENV{TRINITY_HOME}) {
    $ENV{TRINITY_HOME} = "$FindBin::Bin/../../trinityrnaseq/";
}


main: {
    
    my $BASEDIR = cwd();
    
    unless ($ref_trans_fa =~ /^\//) {
        $ref_trans_fa = "$BASEDIR/$ref_trans_fa";
    }
            
            
    unless (-d $OUT_DIR) {
        &process_cmd("mkdir -p $OUT_DIR");
    }
    chdir $OUT_DIR or die "Error, cannot cd to $OUT_DIR";
    

    my $cmd = "";
    
    
    # simulate reads:
    $cmd = "$ENV{TRINITY_HOME}/util/misc/simulate_illuminaPE_from_transcripts.wgsim.pl --transcripts $ref_trans_fa "
        . " --read_length $read_length "
        . " --frag_length $frag_length "
        . " --depth_of_cov 200 "
        . " @ARGV "; # wgsim opts pass-through
        ;
        
    ## todo: add mutation rate info
        
    &process_cmd($cmd);
    
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

