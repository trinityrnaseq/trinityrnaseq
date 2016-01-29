#!/usr/bin/env perl

use strict;
use warnings;

use Cwd;

my $curr_dir = cwd();

foreach my $file (<*_1.*fastq*>, <*_1*.fq*>) {

    $file =~ /^(\S+)_1\.*/ or die "Error, cannot decipher filename $file";
    my $core = $1;

    my $right_fq = $file;
    $right_fq =~ s/_1\./_2\./;
    
    unless (-s $right_fq) {
        die "Error, cannot find Right.fq file corresponding to $file";
    }

    print join("\t", $core, "$curr_dir/$file", "$curr_dir/$right_fq") . "\n";

}


exit(0);


    
