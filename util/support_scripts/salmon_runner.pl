#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;

my $usage = "usage: $0 Trinity.fasta reads.fa\n\n";

my $trin_fa = $ARGV[0] or die $usage;
my $reads_fa = $ARGV[1] or die $usage;





main: {

    my $salmon_index = "$trin_fa.salmon.idx";
    my $cmd = "salmon --no-version-check index -t $trin_fa -i $salmon_index --type quasi -k 25 -p 1";
    &process_cmd($cmd);

    $cmd = "salmon --no-version-check quant -i $salmon_index -l U -r $reads_fa -o salmon_outdir -p 1 --minAssignedFrags 1 ";
    &process_cmd($cmd);
    
    exit(0);

}
