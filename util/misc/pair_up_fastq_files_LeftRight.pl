#!/usr/bin/env perl

use strict;
use warnings;

use Cwd;

my $curr_dir = cwd();

foreach my $file (<*.Left.fq*>, <*.left.fq*>, <*_Left.fq*>, <*_left.fq*>) {

    $file =~ /^(\S+)[\._]Left.fq/i or die "Error, cannot decipher filename $file";
    my $core = $1;

    my $right_fq = $file;
    $right_fq =~ s/([\._])Left/${1}Right/;
    $right_fq =~ s/([\._])left/${1}right/;
    

    unless (-s $right_fq) {
        die "Error, cannot find Right.fq file corresponding to $file";
    }

    print join("\t", $core, "$curr_dir/$file", "$curr_dir/$right_fq") . "\n";

}


exit(0);


    
