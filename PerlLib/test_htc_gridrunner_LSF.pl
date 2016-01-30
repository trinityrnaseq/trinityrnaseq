#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin");

use HTC::GridRunner;

my $config_file = "$FindBin::RealBin/../htc_conf/BroadInst_LSF.test.conf";

main: {

    my @cmds;
    for my $num (1..10) {
        my $cmd = "echo hello $num";
        push (@cmds, $cmd);
    }
    push (@cmds, "this_command_should_fail");
    
    my $grid_runner = new HTC::GridRunner($config_file, "cache_completed_LSF_cmds");

    my $ret = $grid_runner->run_on_grid(@cmds);
    
    exit($ret);
}
    
