#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use FindBin;



my $usage = <<_EOUSAGE_;


##############################################################################################################
#
# Required:
#
#  --sam <string>                   sam-formatted file
#
# Optional:
#
#  --max_insert_size <int>          maximum insert size (default: 500) : used for pair coverage and jaccard
#  --min_insert_size <int>          minimum insert size (default: 100) : used for jaccard
#
#  --jaccard_win_length|-W <int>       default: 100 (requires --jaccard)
#
#  --pseudocounts <int>            default: 1
#  -e                               write extended format, including 'single' and 'both' raw counts
#
##############################################################################################################



_EOUSAGE_

    ;


my $help_flag;

# required
my $sam_file;

# optional
my $max_insert_size = 500;
my $min_insert_size = 100;

my $jaccard_win_length = 100;

my $pseudocounts = 1;
my $extended_flag = 0;


&GetOptions ( 'h' => \$help_flag,

              # required
              'sam=s' => \$sam_file,
              
              # optional
              'max_insert_size=i' => \$max_insert_size,
              'min_insert_size=i' => \$min_insert_size,
              
              'jaccard_win_length|W=i' => \$jaccard_win_length,
              

              'pseudocounts=i' => \$pseudocounts,
              
              'e' => \$extended_flag,

              );


if ($help_flag) {
    die $usage;
}


unless ($sam_file) {
    die $usage;
}

if (@ARGV) {
    die $usage;
}


my $util_dir = "$FindBin::RealBin";

main: {

        
    ## generate the paired coverage info:
    my $cmd = "$util_dir/SAM_to_frag_coords.pl --sam $sam_file "
        . "--max_insert_size $max_insert_size "
        . "--min_insert_size $min_insert_size ";  ## writes file: $sam_file.frag_coords
    
        
    &process_cmd($cmd) unless (-s "$sam_file.frag_coords");
    
    ## compute jaccard coeff info for pair-support
    $cmd = "$util_dir/ordered_fragment_coords_to_jaccard.pl --lend_sorted_frags $sam_file.frag_coords -W $jaccard_win_length "
        . " --pseudocounts $pseudocounts ";
    if ($extended_flag) {
        $cmd .= " -e ";
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
