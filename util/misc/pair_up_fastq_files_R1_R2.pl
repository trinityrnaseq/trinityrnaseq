#!/usr/bin/env perl

use strict;
use warnings;

use Cwd;

my $curr_dir = cwd();

foreach my $file (<*_R1*.fastq*>, <*_R1*.fq*>, <*.R1*.fastq*>, <*.R1*.fq*>) {

    $file =~ /^(\S+)[\._]R1[\._].*f(ast)?q/ or die "Error, cannot decipher filename $file";
    my $core = $1;

    my $right_fq = $file;
    $right_fq =~ s/([\._])R1([\._])/$1R2$2/ or die "Error, cannot convert R1 to R2 in $right_fq";
    
    unless (-s $right_fq) {
        die "Error, cannot find Right.fq file corresponding to $file";
    }

    print join("\t", $core, "$curr_dir/$file", "$curr_dir/$right_fq") . "\n";

}


exit(0);


    
