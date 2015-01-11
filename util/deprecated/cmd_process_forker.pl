#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);
use List::Util qw (shuffle);

my $usage = <<_EOUSAGE_;

#################################################################################
#
#  -c <string>     filename containing a list of commands to execute, one cmd per line.           
#
#  --CPU           Default: 1
#
#  --shuffle       randomly orders the commands before executing them.
#
###################################################################################


_EOUSAGE_

    ;


my $cmds_file = "";
my $CPU = 1;
my $shuffle_flag = 0;

&GetOptions (
             "c=s" => \$cmds_file,
             "CPU=i" => \$CPU,
             "shuffle" => \$shuffle_flag,
             );

unless ($cmds_file) {
    die $usage;
}

my $log_dir = "cmds_log.$$";
mkdir($log_dir) or die "Error, cannot mkdir $log_dir";

## This is very cheap 'parallel-computing' !!!  :)

my $uname = `uname -n`;
chomp $uname;

print "SERVER: $uname, PID: $$\n";

my $ENCOUNTERED_FAILURE = 0;

main: {

    my %job_tracker;
    my @failed_jobs;
    
    my $num_running = 0;
    
    my @unix_cmds;  # just command strings such as you would launch directly from a bash shell.
    my %completed_cmds = ();  # list of completed commands
    my $completed_file_name = $cmds_file.".completed";
    {
        
        ## process list of completed commands
        if(-f ($completed_file_name)){
            open (my $fh, "<", $completed_file_name) or die "Error, cannot open completion file for reading $completed_file_name";
            while (<$fh>) {
                chomp;
                $completed_cmds{$_} = 1;
            }
        }
        ## load all commands in memory... filesystem glitches can otherwise result in premature termination with inadvertent success status.
        open (my $fh, $cmds_file) or die "Error, cannot open file $cmds_file";
        
        while (my $cmd = <$fh>) {
            chomp $cmd;
            if(!$completed_cmds{$cmd}){
                push (@unix_cmds, $cmd);
            }
        }
        close $fh;
    }
    
    if ($shuffle_flag) {
        @unix_cmds = shuffle(@unix_cmds);
    }
    

    my $cmd_counter = 0;
    foreach my $cmd (@unix_cmds) {
        $cmd_counter++;
        
        $num_running++;
        my $child = fork();
        
        if ($child) {
            # parent
            $job_tracker{$cmd_counter} = $cmd;
        }
        else {
            # child:
            my $ret = &run_cmd($cmd_counter, $cmd);
            exit($ret);
        }
        
    
        if ($num_running >= $CPU) {
            wait();
            my $num_finished = &collect_jobs(\%job_tracker, \@failed_jobs, $completed_file_name);
            $num_running -= $num_finished;
            
            ## collect the other finished children to avoid build-up of zombie processes.
            for (1..$num_finished-1) {
                wait();
            }
            
        }
    }
    
    
    ## collect remaining processes.
    while (wait() != -1) { };
    
    &collect_jobs(\%job_tracker, \@failed_jobs, $completed_file_name);
    
    # purge log directory
    `rm -rf $log_dir`;
    
    my $num_failed_jobs = scalar @failed_jobs;
    if (! $ENCOUNTERED_FAILURE) {
        print "\n\nAll $cmd_counter jobs completed successfully! :) \n\n";
        exit(0);
    }
    else {
        unless ($num_failed_jobs == $ENCOUNTERED_FAILURE) {
            print "\n\nError, $ENCOUNTERED_FAILURE jobs failed, but only have recorded $num_failed_jobs ... very strange and unexplained.\n\n";
            # I haven't seen this in my testing, but it's clear that others have, and I haven't figured out why or how yet...  bhaas
        }
        
        # write all failed commands to a file.
        my $failed_cmd_file = "failed_cmds.$$.txt";
        open (my $ofh, ">$failed_cmd_file") or die "Error, cannot write to $failed_cmd_file";
        @failed_jobs = sort {$a->{index}<=>$b->{index}} @failed_jobs;
        foreach my $failed_job (@failed_jobs) {
            print $ofh $failed_job->{cmd} . "\n";
        }
        
        close $ofh;
        
        print "\n\nSorry, $num_failed_jobs of $cmd_counter jobs failed.\n\n"
            . "Failed commands written to file: $failed_cmd_file\n\n";
        exit(1);
    }
}

    
####
sub run_cmd {
    my ($index, $cmd) = @_;

    print "\nRUNNING: $cmd\n";
    
    my $ret = system($cmd);
        
    if ($ret) {
        print STDERR "Error, command: $cmd died with ret $ret";
    }
    
    open (my $log_fh, ">$log_dir/$index.ret") or die "Error, cannot write to log file for $index.ret";
    print $log_fh $ret;
    close $log_fh;


    return($ret);
}


####
sub collect_jobs {
    my ($job_tracker_href, $failed_jobs_aref, $completed_file) = @_;

    my @job_indices = keys %$job_tracker_href;

    my $num_finished = 0;

    foreach my $index (@job_indices) {
        
        my $log_file = "$log_dir/$index.ret";
        
        if (-s $log_file) {
            my $ret_val = `cat $log_file`;
            chomp $ret_val;
            my $job = $job_tracker_href->{$index};
            if ($ret_val == 0) {
                # hurray, job succeded.

                print "SUCCESS[$index]: $job\n";
                # write completed jobs to file
                if(open (my $fh, ">>", $completed_file)){
                    print($fh $job."\n");
                    close($fh);
                }
            }
            else {
                # job failed.
                print "FAILED[$index]: $job\n";
                $ENCOUNTERED_FAILURE++;
                push (@$failed_jobs_aref, {index => $index,
                                           cmd => $job_tracker_href->{$index},
                      });
            }
            
            unlink $log_file;
            $num_finished++;
            delete $job_tracker_href->{$index}; # job done, no need to continue tracking it. (thanks Raj!)
        }
    }

    return($num_finished);
}

