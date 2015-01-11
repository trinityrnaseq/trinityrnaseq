#!/usr/bin/env perl

use strict;
use warnings;

use Cwd;
use File::Basename;

my $usage = "usage: $0 cNumb.graph\n\n";

my $comp = $ARGV[0] or die $usage;

my $TRINITY_HOME = $ENV{TRINITY_HOME} or die "Error, need env var TRINITY_HOME";

main: {


    my $base_comp_name = basename($comp);
    
    my $cwd = cwd();

    unless ($comp =~ /^\//) {
        $comp = "$cwd/$comp";
        
    }

    
    
    
    { ## run old butterfly
        
        my $workdir = "$cwd/__$base_comp_name.oldBfly_dir";

        mkdir($workdir) or die "Error, cannot mkdir $workdir";
        chdir $workdir or die $!;
     
        my $cmd = "ln -s $comp.reads .;  ln -s $comp.out .";
        &process_cmd($cmd);

        $cmd = "java -Xmx4G -jar $TRINITY_HOME/Butterfly/prev_vers/Butterfly_r2013_08_14.jar -N 100000 -L 200 -F 500 -C " . basename($comp) . " --path_reinforcement_distance=75 --max_number_of_paths_per_node=10 -V 15 --stderr 2>&1 | tee log.txt";

        open (my $ofh, ">bfly.cmd") or die $!;
        print $ofh $cmd;
        close $ofh;

        &process_cmd($cmd);

    }
    
    chdir $cwd or die $!;

    
    { ## run new butterfly

        my $workdir = "$cwd/__$base_comp_name.newBfly_dir";
        
        mkdir($workdir) or die $!;
        chdir ($workdir) or die $!;

        my $cmd = "ln -s $comp.reads .; ln -s $comp.out .";
        &process_cmd($cmd);
        

        $cmd = "java -Xmx4G -jar $TRINITY_HOME/Butterfly/Butterfly.jar -N 100000 -L 200 -F 500 -C " . basename($comp) . " --path_reinforcement_distance=75 -V 15 --stderr 2>&1 | tee log.txt";

        open (my $ofh, ">bfly.cmd") or die $!;
        print $ofh $cmd;
        close $ofh;


        &process_cmd($cmd);


    }


    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}


