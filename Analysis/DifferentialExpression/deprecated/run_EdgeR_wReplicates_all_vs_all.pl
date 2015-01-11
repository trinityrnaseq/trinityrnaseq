#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $usage = "usage: $0 samples_described.txt counts.matrix\n\n";

my $samples_descr = $ARGV[0] or die $usage;
my $counts_matrix = $ARGV[1] or die $usage;



=example_samples_descr

Adherent    Adherent_1
Adherent    Adherent_2
Adherent    Adherent_3

Aggregative Aggregative_1
Aggregative Aggregative_2
Aggregative Aggregative_3

Floating    Floating_1
Floating    Floating_2
Floating    Floating_3

=cut


main: {


    my %sample_name_to_reps;
    
    open (my $fh, $samples_descr) or die "Error, cannot open file $samples_descr";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my ($sample_name, $rep_name, @rest) = split(/\s+/);
        push (@{$sample_name_to_reps{$sample_name}}, $rep_name);
    }
    close $fh;

    
    my @samples = keys %sample_name_to_reps;

    for (my $i = 0; $i < $#samples; $i++) {
        my $sample_A = $samples[$i];

        my @reps_A = @{$sample_name_to_reps{$sample_A}};

        for (my $j = $i + 1; $j <= $#samples; $j++) {
            my $sample_B = $samples[$j];
            
            my @reps_B = @{$sample_name_to_reps{$sample_B}};

           
            ## run DESeq
            my $cmd = "$FindBin::Bin/run_EdgeR_wReplicates.pl --counts_matrix $counts_matrix "
                . " --repA_name $sample_A --repA_list \"" . join(",", @reps_A) . "\" "
                . " --repB_name $sample_B --repB_list \"" . join(",", @reps_B) . "\" ";

            &process_cmd($cmd);

        }
    }
     
    exit(0);

}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


