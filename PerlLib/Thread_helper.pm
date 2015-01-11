package Thread_helper;

use strict;
use warnings;
use Carp;
use threads;

=synopsis


    ## here's how you might use it:

    use threads;
    use Thread_helper;

    my $num_simultaneous_threads = 10;
    my $thread_helper = new Thread_helper($num_simultaneous_threads);

    for (1..1000) {

        $thread_helper->wait_for_open_thread();
        
        
        
        my $thread = threads->create(sub{ system "sleep $_";}, int(rand(20)));
        $thread_helper->add_thread($thread);
    }
    $thread_helper->wait_for_all_threads_to_complete();

    my @failures = $thread_helper->get_failed_threads();
    if (@failures) {
        ## examine them...   these are the same threads created above, use use the threads api to access info about them
        ## such as error messages
    }
    else {
        ## all good!
    }


=cut





my $SLEEPTIME = 1;

our $THREAD_MONITORING = 0; # set to 1 to watch thread management


sub new {
    my ($packagename) = shift;
    my ($num_threads) = @_;

    unless ($num_threads && $num_threads =~ /^\d+$/) {
        confess "Error, need number of threads as constructor param";
    }

    my $self = { 
        max_num_threads => $num_threads,
        current_threads => [],
        
        error_threads => [],

        status_success => 0,
        status_error => 0, 

        thread_id_to_command => {},  #tid => cmd 
        thread_id_timing => {}, # tid => { start => val, end => val }
        

    };

    bless ($self, $packagename);

    return($self);
}


sub add_thread {
    my $self = shift;
    my ($thread, $command) = @_;

    my $thread_id = $thread->tid;
   
    $self->{thread_id_to_command}->{$thread_id} = $command;
    $self->{thread_id_timing}->{start}->{$thread_id} = time();
    
    my $num_threads = $self->get_num_threads();
    my $max_num_threads = $self->{max_num_threads};

    if ($num_threads >= $self->{max_num_threads}) {
        
        # if using a wait_for_open_thread() in the client, then this condition won't be met.  Best to keep it in the client rather than here.
        
        print STDERR "- Thread_helper: have $num_threads threads running, waiting for a thread to finish...." if $THREAD_MONITORING;
        $self->wait_for_open_thread();
        print STDERR " done waiting.\n" if $THREAD_MONITORING;
    }
    else {
        print STDERR "- Thread_helper: only $num_threads of max $max_num_threads threads running. Adding another now.\n" if $THREAD_MONITORING;
    }
    
    push (@{$self->{current_threads}}, $thread);

    return;
}


####
sub get_num_threads {
    my $self = shift;
    
    return(scalar @{$self->{current_threads}});
}


####
sub wait_for_open_thread {
    my $self = shift;
    
    if ($self->get_num_threads() >= $self->{max_num_threads}) {
        
        my $waiting_for_thread_to_complete = 1;
        
        my @active_threads;
        
        while ($waiting_for_thread_to_complete) {
            
            @active_threads = ();
            
            my @current_threads = @{$self->{current_threads}};
            foreach my $thread (@current_threads) {
                if ($thread->is_running()) {
                    push (@active_threads, $thread);
                }
                else {
                    $waiting_for_thread_to_complete = 0;
                    $thread->join();
                    my $status;
                    if (my $error = $thread->error()) {
                        my $thread_id = $thread->tid;
                        print STDERR "ERROR, thread $thread_id exited with error $error\n";
                        $self->_add_error_thread($thread);
                        $self->{status_error}++;
                        $status = "ERROR";
                    }
                    else {
                        $self->{status_success}++;
                        $status = "SUCCESS";
                    }
                    my $thread_id = $thread->tid;
                    $self->{thread_id_timing}->{end}->{$thread_id} = time();
                    if ($THREAD_MONITORING) {
                        $self->report_thread_info($thread_id, $status);
                    }
                    
                }
            }
            if ($waiting_for_thread_to_complete) {
                sleep($SLEEPTIME); 
            }
        }
        
        @{$self->{current_threads}} = @active_threads;
        
        
    }

    return;


}



####
sub wait_for_all_threads_to_complete {
    my $self = shift;
    
    print STDERR "\n-now waiting for all threads to complete.\n" if $THREAD_MONITORING;
    
    my @current_threads = @{$self->{current_threads}};
    foreach my $thread (@current_threads) {
        $thread->join();
        my $thread_id = $thread->tid;
        $self->{thread_id_timing}->{end}->{$thread_id} = time();
        my $status = "SUCCESS";
        if (my $error = $thread->error()) {
            print STDERR "ERROR, thread $thread_id exited with error $error\n";
            $self->_add_error_thread($thread);
            $status = "ERROR";
        }
        $self->report_thread_info($thread_id, $status);
    }
    
    @{$self->{current_threads}} = (); # clear them out.
    
    return;

}

####
sub get_failed_threads {
    my $self = shift;

    my @failed_threads = @{$self->{error_threads}};
    return(@failed_threads);
}

####
sub report_thread_info {
    my $self = shift;
    my ($thread_id, $status) = @_;

    my $cmd = $self->{thread_id_to_command}->{$thread_id} || "unknown";
    my $start_time = $self->{thread_id_timing}->{start}->{$thread_id};
    my $end_time = $self->{thread_id_timing}->{end}->{$thread_id};

    print STDERR "Thread($thread_id)\t$status\tCMD: $cmd\tTime to complete: " . ($end_time-$start_time) . " seconds\n";
    
    return;
}


############################
## PRIVATE METHODS #########
############################


####
sub _add_error_thread {
    my $self = shift;
    my ($thread) = @_;

    push (@{$self->{error_threads}}, $thread);

    return;
}


1; #EOM
