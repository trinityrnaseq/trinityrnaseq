package HTC::Base_handler;

use strict;
use warnings;
use Carp;

# basic interface required for implementing a Farm handler.


sub submit_job_to_grid {
    my ($self, $shell_script) = @_;

    confess "abstract method, no implementation here";

    ##########################################################
    # create a command to submit the shell script to the grid.
    ##########################################################



    my $cmd = $self->{config}->{cmd} or confess "Error, need grid cmd template from config";

    #############################################
    # capture and store the corresponding job id
    #############################################
   

    my $job_id_text = `$cmd`;
    # print STDERR "\n$job_id_text\n";
    
    my $ret = $?;
    
    if ($ret) {
        print STDERR "FARMIT failed to accept job: $cmd\n (ret $ret)\n$job_id_text";
        return(-1);
    }
    else {
    
        ## job submitted just fine.
        
        ## get the job ID and log it:
        if ($job_id_text =~ /Job \<(\d+)\>/) {
            my $job_id = $1;
            return($job_id);
            
        }
        else {
            confess "Fatal error, couldn't extract Job ID from submission text: $job_id_text"; 
        }
    }
    
    
}


####
sub job_running_or_pending_on_grid {
    my $self = shift;
    my $job_id = shift;


    confess "Abstract method, example implementation below";
    
    ###########################################################################
    # Given a job_id, query the grid system to determine if it's still running
    ###########################################################################

    unless (defined($job_id)) {
        confess "Error, need job ID as parameter";
    }
    
    # print STDERR "Polling grid to check status of job: $job_id\n";
    
    my $response = `bjobs $job_id`;
    #print STDERR "Response:\n$response\n";

    foreach my $line (split(/\n/, $response)) {
        my @x = split(/\s+/, $line);

        if ($x[0] eq $job_id) {
            my $state = $x[2];
            if ($state eq "DONE" || $state eq "EXIT") {
                return(0);
            }
            else {
                $self->{job_id_to_submission_time}->{$job_id} = time();
                return($state);
            }
        }
    }
    
    print STDERR "-no record of job_id $job_id, setting as state unknown\n";
    return undef; # no status info

}    



1; #EOM
