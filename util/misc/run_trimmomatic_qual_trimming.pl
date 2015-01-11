#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use FindBin;

my $usage = <<__EOUSAGE__;

###############################################################
#
#  --left <string>     left.fq
#  --right <string>    right.fq
#
#    or
#
#  --single <string>   single.fq
#
#  Optional:
#
#  --CPU <int>         default: 4
#
#  --LEADING <int>     default: 15
#  --TRAILING <int>    default: 15
#
#  --MINLEN <int>      default: 36
#
###############################################################

__EOUSAGE__


    ;


my $left;
my $right;
my $single;

my $threads = 4;
my $LEADING = 15;
my $TRAILING = 15;
my $MINLEN = 36;

&GetOptions( 'left=s' => \$left,
             'right=s' => \$right,
             'single=s' => \$single,
             
             'CPU=i' => \$threads,

             'LEADING=i' => \$LEADING,
             'TRAILING=i' => \$TRAILING,
             'MINLEN=i' => \$MINLEN,
             );


=trimmomatic

java -jar /seq/regev_genome_portal/SOFTWARE/BIN/trimmomatic.jar PE -threads {__THREADS__} -phred33 \
{__LEFT_FQ__} {__RIGHT_FQ__} \
{__LEFT_FQ__}.P.qtrim.fq {__LEFT_FQ__}.U.qtrim.fq  \
{__RIGHT_FQ__}.P.qtrim.fq {__RIGHT_FQ__}.U.qtrim.fq \
  LEADING:15 TRAILING:15 MINLEN:36 2> {__LOCAL_ANALYSIS_DIR__}/trimmomatic.log.stats

=cut

    ;

unless ( ($left && $right) || $single) {
    die $usage;
}

main: {

    my $cmd;

    if ($left && $right) {
    
        $cmd = "java -jar $FindBin::Bin/../trinity-plugins/Trimmomatic/trimmomatic.jar PE -threads $threads -phred33 "
            . " $left $right "
            . " $left.P.qtrim.fq $left.U.qtrim.fq "
            . " $right.P.qtrim.fq $right.U.qtrim.fq "
            . " LEADING:$LEADING TRAILING:$TRAILING MINLEN:$MINLEN ";
    }
    else {
        
        $cmd = "java -jar $FindBin::Bin/../trinity-plugins/Trimmomatic/trimmomatic.jar SE -threads $threads -phred33 "
            . " $single "
            . " $single.qtrim.fq "
            . " LEADING:$LEADING TRAILING:$TRAILING MINLEN:$MINLEN ";
        
    }

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
