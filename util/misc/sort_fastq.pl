#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fastq_reader;

my $usage = "usage: $0 file.fastq\n\n";

my $fastq_file = $ARGV[0] or die $usage;

main: {

    my $fastq_reader = new Fastq_reader($fastq_file);

    my $tab_tmp_file = "$fastq_file.tab";

    open (my $ofh, ">$tab_tmp_file") or die $!;

    while (my $fq_entry = $fastq_reader->next()) {

        my $fq_record = $fq_entry->get_fastq_record();
        chomp $fq_record;
        
        my @lines = split(/\n/, $fq_record);
        
        print $ofh join("\t", @lines) . "\n";

    }
    close $ofh;

    ## sort by read name

    my $cmd = "sort -S 4G -T . -k1,1 $tab_tmp_file > $tab_tmp_file.sort";
    &process_cmd($cmd);
    
    unlink($tab_tmp_file);
    
    ## convert back to fastq file format

    $cmd = "cat $tab_tmp_file.sort | sed s/\\\\t/\\\\n/g > $fastq_file.sorted.fq";
    &process_cmd($cmd);

    unlink("$tab_tmp_file.sort");
    
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

