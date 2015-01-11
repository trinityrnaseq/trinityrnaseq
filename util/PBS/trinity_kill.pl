#!/usr/bin/env perl

use strict;
use warnings;

my $me = $ENV{'USER'};
my @jobs_running_ln = `qstat -u $ENV{'USER'}`;
if (!@jobs_running_ln || scalar(@jobs_running_ln)<1){
	print "No running jobs for user $me\n";
	exit();
}

my %jobs_running;
foreach my $job_ln (@jobs_running_ln){
	$job_ln=~/^(\S+)/;
	$jobs_running{$1}=1 if $1;
}

my $dir=$ARGV[0] ? $ARGV[0] : '.';
if ($dir && -d $dir){
	# read jobnumbers in current directory (and directory passed as variable) and kill them (start at the last job and move to oldest)
	if (-s $dir."/jobnumbers.out"){
		my @job_sub = `tac $dir/jobnumbers.out`;
		chomp(@job_sub);
		foreach my $job (@job_sub){
			system("qdel -W force $job") if $jobs_running{$job};
			# twice to make sure
			system("qdel -W force $job >/dev/null 2>/dev/null") if $jobs_running{$job};
		}
		unlink($dir."/jobnumbers.out");
	}
	else{
		print "No previous jobs found\n";
	}
}

