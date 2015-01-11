package HTC::GridRunner;

use strict;
use warnings;
use Carp;
use List::Util qw (shuffle);
use FindBin;
use Cwd;

use HTC::FarmIt;
use HTC::LSF_handler;
use HTC::SGE_handler;
use HTC::SLURM_handler;
use HTC::PBS_handler;

BEGIN {
    unless ($ENV{TRINITY_HOME}) {
        confess "ERROR, must have env var TRINITY_HOME set to the base installation directory of Trinity ";
    }

    unless ($ENV{HOSTNAME}) {
        if ($ENV{HOST}) {
            $ENV{HOSTNAME} = $ENV{HOST};
        }
        else {
            $ENV{HOSTNAME} = `hostname`;
        }
    }
    
}

####
sub new {
    my $packagename = shift;
    my $config_file = shift;
    my $cache_completed_cmds_file = shift;

    unless ($config_file) {
        confess "Error, need config file as parameter";
    }
    
    my %config;
    open (my $fh, $config_file) or confess "Error, cannot open file: $config_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        
        my ($tok, $val) = split(/=/, $_, 2);
        $config{$tok} = $val;
    }
    close $fh;

    use Data::Dumper;
    print STDERR Dumper(\%config);

    my $grid_type = $config{grid} or confess "Error, grid type not specified in config file: $config_file";
    unless ($grid_type =~ /^(LSF|SGE|SLURM|PBS)$/) {
        confess "Error, grid type: $grid_type is not currently supported";
    }
    
    my $handler;
    if ($grid_type eq "LSF") {
        $handler = HTC::LSF_handler->new(\%config);
    }
    elsif ($grid_type eq "SGE") {
        $handler = HTC::SGE_handler->new(\%config);
    }
    elsif ($grid_type eq "SLURM") {
        $handler = HTC::SLURM_handler->new(\%config);
    }
    elsif ($grid_type eq "PBS") {
        $handler = HTC::PBS_handler->new(\%config);
    }
    else {
        confess "Error, grid type: $grid_type is not supported";
    }
    

    unless ($cache_completed_cmds_file) {
        $cache_completed_cmds_file = "$$.completed_cmds.cache";
    }

    my $self = { handler => $handler,
                 config => \%config,
                 cache_completed_cmds_file => $cache_completed_cmds_file,
             };
    
    bless($self, $packagename);

    return($self);
}
    


####
sub run_on_grid {
    my $self = shift;
    my @cmds = @_;
    
    # check for earlier-completed commands from a previous run. Don't rerun them.
    my $cache_completed_cmds_file = $self->{cache_completed_cmds_file};
    my %completed_cmds;
    if ($cache_completed_cmds_file && -s $cache_completed_cmds_file) {
        open (my $fh, $cache_completed_cmds_file) or confess "Error, cannot open file: $cache_completed_cmds_file";
        while (<$fh>) {
            my $cmd = $_;
            chomp $cmd;
            $completed_cmds{$cmd} = 1;
        }
        close $fh;
        
        my @unfinished_cmds;
        my $count_finished_before = 0;
        foreach my $cmd (@cmds) {
            if (! exists $completed_cmds{$cmd}) {
                push (@unfinished_cmds, $cmd);
            }
            else {
                $count_finished_before++;
            }
        }
        if ($count_finished_before) {
            print STDERR "-note, $count_finished_before commands already completed successfully. Skipping them here.\n";
            print STDERR "-there are " . scalar(@unfinished_cmds) . " cmds left to run here.\n";
            @cmds = @unfinished_cmds;
        }
    }
    
    
    @cmds = shuffle @cmds;

    my %params = ( cmds => \@cmds,
                   handler => $self->{handler},
                   );

    

    if ($cache_completed_cmds_file) {
        open (my $ofh, ">>$cache_completed_cmds_file") or die "Error, cannot append to file $cache_completed_cmds_file";
        $params{cache_success_cmds_ofh} = $ofh;
    }
    
    if (my $max_nodes = $self->{config}->{max_nodes}) {
        $params{max_nodes} = $max_nodes;
    }
    if (my $cmds_per_node = $self->{config}->{cmds_per_node}) {
        $params{cmds_per_node} = $cmds_per_node;
    }

    my $farmer = new HTC::FarmIt(\%params);
    
    $farmer->submit_jobs();
    
    
    my @failed_cmds = $farmer->get_failed_cmds();

    
    if (@failed_cmds) {
        my $num_failed_cmds = scalar(@failed_cmds);
        
        my @cmds_remaining;
        foreach my $failed_cmd_info (@failed_cmds) {
            push (@cmds_remaining, $failed_cmd_info->{cmd});
        }

        print STDERR "$num_failed_cmds commands failed during grid computing.\n";
        print STDERR "Failed cmds:\n" . join("\n", @cmds_remaining) . "\n";
        

        ## try running them via parafly
        print STDERR "\n\nTrying to run them using parafly...\n\n";
        return($self->run_parafly(@cmds_remaining));
        
    }
    else {
        print "All commands completed successfully on the computing grid.\n";
        return(0);
    }
}

####
sub run_parafly {
    my $self = shift;
    my (@cmds) = @_;

    my $cache_file = $self->{cache_completed_cmds_file};

    my $cmds_file = "$cache_file.failed_for_parafly.txt";
    open (my $ofh, ">$cmds_file") or die "Error, cannot write to file $cmds_file";
    foreach my $cmd (@cmds) {
        print $ofh $cmd . "\n";
    }

    my $num_cpus = $ENV{OMP_NUM_THREADS} || 10;

    my $cmd;
    if (my $trin_home = $ENV{TRINITY_HOME}) {
        $cmd = "$trin_home/trinity-plugins/parafly/bin/ParaFly";
    }
    else {
        $cmd = "ParaFly";
    }

    $cmd .= " -c $cmds_file -CPU $num_cpus -v -shuffle -failed_cmds $cmds_file.FAILED_DURING_PARAFLY";
    
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret: $ret";
    }
    
    return(0);
}
    
1;

    

