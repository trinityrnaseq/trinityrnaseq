#!/usr/bin/env perl

package main;
our $SEE;

package HPC::FarmIt;

use strict;
use warnings;
use Cwd;
use File::Basename;
use Carp;
use Sys::Hostname;

my $WAITTIME = 15;
my $RETVAL_BIN_SIZE = 1000;
#my $RESORT_TO_POLLING_TIME = 60*60; # 1 hour (60 min * 60 seconds)
my $RESORT_TO_POLLING_TIME = 15*60; # 15 minutes

my %JOB_ID_TO_PREVTIME;

my $DEFAULT_MAX_NODES = 500;
my $DEFAULT_CMDS_PER_NODE = 10;
my $TMPDIR = $ENV{TMPDIR} || "/local/scratch";

my $MADE_LOGDIR_HERE_FLAG = 0;

=description

    Here\'s some sample code that demonstrates typical usage;


my @cmds = (list of unix commands to execute)


use List::Util qw (shuffle);
@cmds = shuffle (@cmds);

   
my $farmer = new Farmit({cmds=>\@cmds,
                         handler => HPC::LSF->new(...), or HPC::SGE->new(...)

                             # optional, have defaults:
                             log_base_dir => cwd(),
                             cmds_per_node => 10,
                             max_nodes => 500,
                             
                            
                            }
                       );

$farmer->submit_jobs();

my $total_cmds = scalar (@cmds);

if (my @failed_cmds = $farmer->get_failed_cmds()) {
    
    my $num_failed_cmds = scalar (@failed_cmds);
    print "Sorry, $num_failed_cmds of $total_cmds failed = " . ($num_failed_cmds
 / $total_cmds * 100) . " % failure.\n";

    open (my $failed_fh, ">failed_cmds") or die $!;
    foreach my $failed_cmd (@failed_cmds) {
        my $cmd = $failed_cmd->{cmd};
        my $ret = $failed_cmd->{ret};
        my $index = $failed_cmd->{index};
        print $failed_fh "# cmd($index)\n$cmd\nret($ret)\n\n";
    }
    close $failed_fh;
}
else {
    print "All $total_cmds completed successfully.\n\n";
}

$farmer->clean_logs();


=cut
    
    ;

sub new {
    my $packagename = shift;
    
    my $params_href = shift;

    umask(0000);
    
    my $host = hostname
    # required:
    my $cmds_list_aref = $params_href->{cmds} or confess "No commands specified";
    my $handler = $params_href->{handler} or confess "Need handler specified";
    

    # optional, have defaults
    my $max_nodes = $params_href->{max_nodes} || $DEFAULT_MAX_NODES;
    my $cmds_per_node = $params_href->{cmds_per_node} || $DEFAULT_CMDS_PER_NODE;
    
    my $log_base_dir = $params_href->{log_base_dir} || cwd();
    
    my $log_dir .= "$log_base_dir/farmit.J$$.$host.$$." . time();
    unless (-d $log_dir) {
        my $ret = system("mkdir -p $log_dir");
        if ($ret) { die "Error, cannot mkdir $log_dir";}
        $MADE_LOGDIR_HERE_FLAG = 1;
    }
        
    ## write commands listing:
    open (my $fh, ">$log_dir/cmds_list.txt") or die $!;
    my $index = 0;
    foreach my $cmd (@$cmds_list_aref) {
        print $fh "index($index)\t$cmd\n";
        $index++;
    }
    close $fh;
    

    
    # finish logdir setup and object creation.
    my $cmds_dir = "$log_dir/cmds";
    my $retvals_dir = "$log_dir/retvals";
    my $monitor_dir = "$log_dir/monitor";
    foreach my $dir ($cmds_dir, $retvals_dir, $monitor_dir) {
        mkdir $dir or die "Error, cannot mkdir $dir";
    }
    
    my $num_cmds = scalar (@$cmds_list_aref);
    my $self = {
        num_cmds => $num_cmds,
        cmds_list => $cmds_list_aref,
        log_dir => $log_dir,
        cmds_dir => $cmds_dir,
        retvals_dir => $retvals_dir,
        monitor_dir => $monitor_dir,
        max_nodes => $max_nodes,
        cmds_per_node => $cmds_per_node,
        retvalues => [], 
        nodes_in_progress => {},
        handler => $handler,
        
        
        
                
    };
    
    if (my $ofh = $params_href->{cache_success_cmds_ofh}) {
        my $curr_fh = select($ofh);
        $| = 1; # flush always
        select $curr_fh;
        $self->{cache_success_cmds_ofh} = $ofh; 
    }
    
    
    bless ($self, $packagename);
    return ($self);
}

####
sub submit_jobs {
    my $self = shift;
   
    $self->_write_pid_file();
    
    
    my $max_nodes = $self->{max_nodes};
    my $num_cmds = $self->{num_cmds};

    
    my $num_cmds_launched = 0;
    my $num_nodes_used = 0;
    

    while ($num_cmds_launched < $num_cmds) {
        $num_cmds_launched = $self->_submit_job($num_cmds_launched);
        $num_nodes_used = $self->_get_num_nodes_used();
        print STDERR "\r  CMDS: $num_cmds_launched / $num_cmds  [$num_nodes_used/$max_nodes nodes in use]   ";
        if ($num_nodes_used >= $max_nodes) {
            my $num_nodes_finished = $self->_wait_for_completions();
            $num_nodes_used -= $num_nodes_finished;
        }
    }
    
    print STDERR "\n* All cmds submitted to grid.  Now waiting for them to finish.\n";
    ## wait for rest to finish
    while (my $num_nodes_finished = $self->_wait_for_completions()) { 
        $num_nodes_used -= $num_nodes_finished;
        print STDERR "\r  CMDS: $num_cmds_launched / $num_cmds  [$num_nodes_used/$max_nodes nodes in use]   ";
    };
    
    print STDERR "\n* All nodes completed.  Now auditing job completion status values\n";
    

    my $retvals_aref = $self->{retvalues};
    
    my $num_successes = 0;
    my $num_failures = 0;
    my $num_unknown = 0;

    foreach my $retval (@$retvals_aref) {
        if ($retval =~ /\d+/) {
            if ($retval == 0) {
                $num_successes++;
            } else {
                $num_failures++;
            }
        } else {
            $num_unknown++;
        }
    }


    $self->_write_result_summary($num_successes, $num_failures, $num_unknown);
    
    if ($num_successes == $num_cmds) {
        print "All commands completed successfully.\n";
    } else {
        print "Failures encountered:\n"
            . "num_success: $num_successes\tnum_fail: $num_failures\tnum_unknown: $num_unknown\n";
    }




    print "Finished.\n\n";
}


####
sub _submit_job {
    my $self = shift;
    my $num_cmds_launched = shift;
    
    my $num_cmds = $self->{num_cmds};
    my $log_dir = $self->{log_dir};
    my $retvals_dir = $self->{retvals_dir};
    my $cmds_dir = $self->{cmds_dir};
    my $cmds_per_node = $self->{cmds_per_node};
    my $cmds_list_aref = $self->{cmds_list};
    my $monitor_dir = $self->{monitor_dir};

    my $orig_num_cmds_launched = $num_cmds_launched;
    
    
    my $shell_script = "$cmds_dir/J$$.S${num_cmds_launched}.sh";
    open (my $fh, ">$shell_script") or die $!;
    print $fh "#!/bin/sh\n\n";
    
    &_write_minimal_environment($fh);
    
    my $num_cmds_written = 0;

    my $monitor_started = "$monitor_dir/$num_cmds_launched.started";
    my $monitor_finished = "$monitor_dir/$num_cmds_launched.finished";

    my @cmd_indices_prepped;
    
    while ($num_cmds_launched < $num_cmds && $num_cmds_written < $cmds_per_node) {
        my $next_cmd_index = $num_cmds_launched; #always one less than the current index
        my $cmd_string = $cmds_list_aref->[ $next_cmd_index ];
        
        push (@cmd_indices_prepped, $next_cmd_index);
        
		my $retval_bin = int($next_cmd_index / $RETVAL_BIN_SIZE);
        
		my $retval_subdir = "$retvals_dir/$retval_bin";
		unless (-d $retval_subdir) {
			mkdir $retval_subdir or die "Error, cannot mkdir $retval_subdir";
		}

        print $fh "## Command index $next_cmd_index\n"
            . "touch $monitor_started\n"
            . "$cmd_string\n"
            . 'echo $? >> ' . "$retval_subdir/entry_$next_cmd_index.ret\n\n";
        
        $num_cmds_launched++;
        $num_cmds_written++;
    }
    
    print $fh "\n" 
        . "rm -f $monitor_started\n"
        . "touch $monitor_finished\n"
        . "\n" 
        . "exit 0\n\n";
    
    
    close $fh;
    chmod (0775, $shell_script);
    
    print "Submitting: $shell_script to farmit\n" if $SEE;
    
    my $job_id = $self->{handler}->submit_job_to_grid($shell_script);
    
    if ($job_id < 0) {
        print STDERR "FARMIT failed to accept job.  Will try again shortly.\n";
        
        unlink $shell_script; # cleanup, try again later
        
        sleep(2*60); # sleep 2 minutes for now.  Give the system time to recuperate if a problem exists
        return ($orig_num_cmds_launched);
        
    }
    else {

        # submitted just fine. 
        my $orig_shell_script_name = $shell_script;
        
        $shell_script = basename($shell_script);
        open (my $logdir_jobsfh, ">>$log_dir/job_ids.txt") or die "Error, cannot open file $log_dir/job_ids.txt";
        ## get the job ID and log it:
        print $logdir_jobsfh "$job_id\t$shell_script\n";
        my $monitor_href = $self->{nodes_in_progress};
        $monitor_href->{$monitor_finished} = $job_id;
            
        $self->{job_id_to_cmd_indices}->{$job_id} = \@cmd_indices_prepped;
        $self->{job_id_to_submission_time}->{$job_id} = time();
        $self->{job_id_to_shell_script}->{$job_id} = $orig_shell_script_name;
        

        close $logdir_jobsfh;
        
        return ($num_cmds_launched);
    }
    
}

sub _get_num_nodes_used {
    my $self = shift;
    my $num_nodes_used = scalar (keys %{$self->{nodes_in_progress}});
    
    print "Num nodes currently in use: $num_nodes_used\n" if $SEE;
    
    return ($num_nodes_used);
}



sub get_failed_cmds {
    my $self = shift;
    my $retvalues_aref = $self->{retvalues};
    my $cmds_list_aref = $self->{cmds_list};

    my @failed_cmds;
    for (my $i = 0; $i <= $#$retvalues_aref; $i++) {
        my $retval = $retvalues_aref->[$i];
        if ($retval) {
            push (@failed_cmds, 
                  { cmd => $cmds_list_aref->[$i],
                    ret => $retval,
                } );
        }
    }
    return (@failed_cmds);
}





sub _wait_for_completions {
    my $self = shift;
    
    print "sub _wait_for_completions()\n" if $SEE;
    
    my $nodes_in_progress_href = $self->{nodes_in_progress};
    
    my $seen_finished = 0;

    my @done;
    while (! $seen_finished) {
                
        ## check to see if there are any jobs remaining:
        if ($self->_get_num_nodes_used() == 0) {
            ## no jobs in the queue
            print "no nodes in use; exiting wait.\n" if $SEE;
            return (0);
        }
        
        ## check for finished jobs
        foreach my $monitor_file (keys %$nodes_in_progress_href) {
            if (-e $monitor_file) {
                push (@done, $monitor_file);
                $seen_finished = 1;
            }
            else {
                ## try polling the grid directly based on the job id
                my $job_id = $nodes_in_progress_href->{$monitor_file};
                
                my $time_launched = $self->{job_id_to_submission_time}->{$job_id};
                my $current_time = time();
                
                ## see if an hour has passed
                if ($current_time - $time_launched >= $RESORT_TO_POLLING_TIME) {
                    ## poll the system directly:
                    if (! $self->job_running_or_pending_on_grid($job_id)) {
                        
                        push (@done, $monitor_file);
                        $seen_finished = 1;
                        
                    }
                    else {
                        ## reset submission time to delay next polling time
                        $self->{job_id_to_submission_time}->{$job_id} = time();
                    }
                    
                }
            }
            
        }
        if ($seen_finished) {
            foreach my $monitor_file (@done) {
                my $job_id = $nodes_in_progress_href->{$monitor_file};
                print "job[$job_id]: $monitor_file is finished.\n" if $SEE;
                
                my $shell_script = $self->{job_id_to_shell_script}->{$job_id}; 
                
                my @cmd_indices = @{$self->{job_id_to_cmd_indices}->{$job_id}};
                my $all_OK = 1;
                foreach my $cmd_index (@cmd_indices) {
                    my $retval_file = $self->_get_ret_filename($cmd_index);
                    if (-s $retval_file) {
                        open (my $fh, $retval_file) or die $!;
                        my $retval_string = <$fh>;
                        $retval_string =~ s/\s//g;
                        $self->{retvalues}->[$cmd_index] = $retval_string;
                        close $fh;
                        # write success if success
                        if ($retval_string == 0) {
                            # command succeeded.
                            if (my $ofh = $self->{cache_success_cmds_ofh}) {
                                my $cmd = $self->{cmds_list}->[$cmd_index];
                                print $ofh "$cmd\n";
                            }
                        }
                        else {
                            # at least one failure
                            $all_OK = 0;
                        }
                        unlink($retval_file); # house cleaning
                    } else {
                        $self->{retvalues}->[$cmd_index] = "FILE_NOT_EXISTS";
                    }
                }
                
                delete $nodes_in_progress_href->{$monitor_file}; #remove from queue
                delete $self->{job_id_to_cmd_indices}->{$job_id};
                delete $self->{job_id_to_submission_time}->{$job_id};
            
                unlink($monitor_file); # house cleaning
                
                if ($all_OK) {
                    # lets not retain the shell script and stderr and stdout.  Only retain them for failures.
                    unlink($shell_script);
                    unlink("$shell_script.stderr", "$shell_script.stdout");
                    delete $self->{job_id_to_shell_script}->{$job_id};
                    
                    print "OK.  Trying to delete $shell_script and .stderr and .stdout.\n" if $SEE;
                    
                }
                else {
                    print "Not OK!!!\n" if $SEE;
                }
                
            }
            return (scalar (@done)); #num jobs completed
        } 
        else {
            ## wait a while and check again
            print "waiting for jobs to finish.\n" if $SEE;
            sleep($WAITTIME);
        }
    }
}


sub _write_pid_file {
    my $self = shift;
    my $log_dir = $self->{log_dir};
    my $host = hostname;
    open (my $fh, ">$log_dir/$host.pid") or die $!;
    print $fh $$;
    close $fh;
}


sub _write_result_summary {
    my ($self, $num_successes, $num_failures, $num_unknown) = @_;
    my $status = ($num_failures == 0 && $num_unknown == 0) ? "success" : "failure"; 
   
    $self->{status} = $status;
    $self->{num_failures} = $num_failures;
    $self->{num_successes} = $num_successes;
    $self->{num_unknown} = $num_unknown;
    
    my $log_dir = $self->{log_dir};
    open (my $fh, ">$log_dir/farmit.finished.$status") or die $!;
    print $fh "num_successes: $num_successes\n"
        . "num_failures: $num_failures\n"
        . "num_unknown: $num_unknown\n";
    close $fh;
    
}

sub clean_logs {
    my $self = shift;
    my $log_dir = $self->{log_dir};
    
    if ($MADE_LOGDIR_HERE_FLAG) {
        my $cmd = "rm -rf $log_dir";
        system $cmd;
        return ($?);
    }
    else {
        print STDERR "** WARNING ** : not removing log_dir: $log_dir because not created by this process.  Security issue...";
        return(-1);
    }
}


sub _write_minimal_environment {
    my ($fh) = @_;

    print $fh <<_EOFENV_;

## add any special environment settings

echo HOST: \$HOSTNAME
echo HOST: \$HOSTNAME >&2

_EOFENV_

;

    return;
    
}



####
sub _get_ret_filename {
    my $self = shift;
    my ($cmd_index) = @_;

    my $retvals_dir = $self->{retvals_dir};

    my $retval_bin = int ($cmd_index / $RETVAL_BIN_SIZE);
    my $retval_file = $retvals_dir . "/$retval_bin/entry_$cmd_index.ret";
    
    return($retval_file);
}



####
sub job_running_or_pending_on_grid {
    my $self = shift;
    my ($job_id) = @_;
    
    if (time() - $self->{job_id_to_submission_time}->{$job_id} < $RESORT_TO_POLLING_TIME) {
        return("TOO_SOON");
    }

    # have handler check using handler-specific method.
    return($self->{handler}->job_running_or_pending_on_grid($job_id));
    
}


1; #EOM




