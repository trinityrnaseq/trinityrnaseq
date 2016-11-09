#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $help_flag;

my $CPU = 4;
my $output_dir = "DetonateData";

my $species_opts = "mouse|human|fission_yeast";

my $usage = <<__EOUSAGE__;

#######################################################################################
#
#  (note, must set env var DETONATE_HOME to installation directory of DETONATE software)
#
#  --reads <string>       if paired-end, list as comma-delimited:  "left.fq,right.fq"
#
#  --target <string>      target fasta file (eg. Trinity.fasta)
#
#  --frag_len <string>    frag length. If SE data, it's the read length.
#
#  --species <string>     $species_opts  (param files provided at: rsem-eval/true_transcript_length_distribution/)
#
#  optional:
#
#  --SS_lib_type <string>  R, F, RF, or FR
#
#  --threads <int>        number of threads to use in multi-threading (default: $CPU)
#
#  --output_dir <string>  default: $output_dir
#
########################################################################################


__EOUSAGE__

    ;


my $reads;
my $target;
my $frag_len;
my $SS_lib_type;
my $species;

&GetOptions ( 'h' => \$help_flag,
              
              'reads=s' => \$reads,
              'target=s' => \$target,
              'frag_len=i' => \$frag_len,
              'SS_lib_type=s' => \$SS_lib_type,
              'threads=i' => \$CPU,
              'output_dir=s' => \$output_dir,
              'species=s' => \$species,
    );


main: {

    my $detonate_home_dir = $ENV{DETONATE_HOME} or die "Error, must set env var DETONATE_HOME to its installation directory";

    unless ($species =~ /$species_opts/) {
        die "Error, species $species not supported, only $species_opts";
    }

    my $cmd = "$detonate_home_dir/rsem-eval-calculate-score -p $CPU "
     . " --transcript-length-parameters $detonate_home_dir/rsem-eval/true_transcript_length_distribution/$species.txt "
     . " $reads $target $output_dir $frag_len ";

    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    else {
        print STDERR "Done.\n";
    }
    exit(0);
    
}
