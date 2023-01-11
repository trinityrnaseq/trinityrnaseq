#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 samples_file dir/containing/bams\n\n";

my $samples_file = $ARGV[0] or die $usage;
my $bam_dir_path = $ARGV[1] or die $usage;


my %tissue_id_to_bams;

open(my $fh, $samples_file) or die $!;
while(<$fh>) {
    chomp;
    my ($tissue_type, $sample_id, @rest) = split(/\t/);
    my @bams = <$bam_dir_path/$sample_id.*.bam>;
    unless (scalar(@bams) == 1) {
        die "Error, cannot find bam corresping to $bam_dir_path/$sample_id.*.bam ";
    }
    push (@{$tissue_id_to_bams{$tissue_type}}, $bams[0]);
}
close $fh;

foreach my $tissue_id (keys %tissue_id_to_bams) {
    
    my @bams = @{$tissue_id_to_bams{$tissue_id}};
    
    my $cmd = "samtools merge $tissue_id.bam @bams";
    print "$cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, $cmd died with ret $ret";
    }
}

exit(0)


