#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use HTC::GridRunner;
use List::Util qw (shuffle);

use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<_EOUSAGE_;

################################################################
# Required:
#
#  -c <string>        file containing list of commands
#  --grid_conf|G <string>   grid config file
#
####################################################################

_EOUSAGE_

	;


my $grid_conf_file;
my $cmd_file;
my $help_flag;

&GetOptions ( 'h' => \$help_flag,
			  'c=s' => \$cmd_file,
		
              'grid_conf|G=s' => \$grid_conf_file,
             
              );


unless ($cmd_file && $grid_conf_file) { 
	die $usage;
}

if ($help_flag) {
	die $usage;
}



## add Parafly to path
$ENV{PATH} .= ":" . "$FindBin::Bin/../trinity-plugins/parafly/bin/";


main: {

	my $uname = `uname -n`;
	chomp $uname;

	print "SERVER: $uname, PID: $$\n";
	
    
    open (my $fh, $cmd_file) or die "Error, cannot open $cmd_file";
    my @cmds;

    while (<$fh>) {
        chomp;
        if (/\w/) {
            push (@cmds, $_);
        }
    }
    close $fh;

    @cmds = shuffle @cmds;  ## to even out load on grid nodes.  Some may topload their jobs!

    my $cache_file = "$cmd_file.htc-cache_success";
    
    my $grid_runner = new HTC::GridRunner($grid_conf_file, $cache_file);
    my $ret = $grid_runner->run_on_grid(@cmds);
        
    if ($ret) {
        
        print STDERR "Error, not all commands could complete successfully... cannot continue.";
        
        exit(1);
    }
    else {
        ## all good
        exit(0);
    }
}


    
