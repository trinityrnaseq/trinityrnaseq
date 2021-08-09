#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Process_cmd;

my $usage = "\n\n\tusage: $0 Trinity.fasta reads.fa [threads=1]\n\n";

my $trin_fa = $ARGV[0] or die $usage;
my $reads_fa = $ARGV[1] or die $usage;
my $CPU = $ARGV[2] || 1;




main: {

    my $salmon_index = "$trin_fa.salmon.idx";
    my $salmon_stderr = "_salmon.$$.stderr";
    my $cmd = "salmon --no-version-check index -t $trin_fa -i $salmon_index -k 25 -p $CPU  > $salmon_stderr 2>&1";
    &run_cmd_capture_stderr($cmd, $salmon_stderr);

    
    $cmd = "salmon --no-version-check quant -i $salmon_index -l U -r $reads_fa -o salmon_outdir -p $CPU --minAssignedFrags 1 --validateMappings > $salmon_stderr 2>&1";
    &run_cmd_capture_stderr($cmd, $salmon_stderr);

    unlink($salmon_stderr);
    
    exit(0);

}


####
sub run_cmd_capture_stderr {
    my ($cmd, $stderr_file) = @_; 
    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        my $errmsg = `cat $stderr_file`;
        die "Error, cmd: $cmd failed with msg: $errmsg $@";
    }
    return;
}
