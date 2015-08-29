package HPC::PBS_handler;

use strict;
use warnings;
use base qw(HPC::Base_handler);
use Carp;
use Cwd;

####
sub new {
    my $packagename = shift;
    my $config_href = shift;
    
    unless (ref $config_href) {
        confess "Error, need config href as param";
    }
    
    
    my $self = {config => $config_href};
    
    bless($self, $packagename);

    return($self);
}

####
sub submit_job_to_grid {
    my $self = shift;
    my $shell_script = shift;

    unless (-s $shell_script) {
        confess "Error, need shell script to submit as parameter.";
    }

    ## submit the command, do any additional job administration as required, such as capturing job ID

    my $cmd = $self->{config}->{cmd} or confess "Error, need cmd from config file";
    
    $cmd .= " -d .  -j oe -N $shell_script ";

    $cmd .= " $shell_script 2>&1 ";
    

    print STDERR "CMD: $cmd\n";

    my $job_id_text = `$cmd`;
    
    
    
    
    my $ret = $?;
    
    print STDERR "\nPBS job (ret=$ret), text:: $job_id_text\n";
    

    if ($ret) {
        print STDERR "FARMIT failed to accept job: $cmd\n (ret $ret)\n$job_id_text\n";
        return(-1);
    }
    else {
    
        ## job submitted just fine.
        
        ## get the job ID and log it:
        if ($job_id_text =~ /(\d+)(\.[0-9a-zA-Z])*/) {
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

    unless (defined($job_id)) {
        confess "Error, need job ID as parameter";
    }
    
    # print STDERR "Polling grid to check status of job: $job_id\n";
    my $response = `qstat`;
    #print STDERR "Response:\n$response\n";

    foreach my $line (split(/\n/, $response)) {
        my @x = split(/\s+/, $line);

        if ($x[0] eq $job_id) {
            my $state = $x[4];
            
            $self->{job_id_to_submission_time}->{$job_id} = time();
            return($state);
            
        }
    }
    
    print STDERR "-no record of job_id $job_id, setting as state unknown\n";
    return undef; # no status info

}    


1; #EOM
