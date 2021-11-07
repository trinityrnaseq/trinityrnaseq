#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $usage = "\n\nusage: $0 --prefix \$prefix fileA.sra [fileB.sra ...]\n\n\n";


my $prefix;

&GetOptions ('prefix=s' => \$prefix);

my @sra_files = @ARGV;

unless (@ARGV) {
    die $usage;
}

unless (@sra_files) {
    die $usage;
}

my $fastq_dump_path = `sh -c "command -v fastq-dump"`;
unless ($fastq_dump_path && $fastq_dump_path =~ /\w/) {
    die "Error, cannot find 'fastq-dump' utility in your PATH. Be sure you have SRA toolkit installed and fastq-dump in your PATH setting. ";
}

my @core_names;
foreach my $sra_file (@sra_files) {

    my $core_name = $sra_file;
    $core_name =~ s/\.sra$//;
    $core_name = basename($core_name);
    push (@core_names, $core_name);
    
    my $cmd = "fastq-dump --defline-seq '@\$sn[_\$rn]/\$ri' --split-files $sra_file";
    &process_cmd($cmd) unless (-s "${core_name}_1.fastq");
    
}

if ($prefix) {
    
    my @final_cmds;
    my @tmp_files;
    for my $end ("1", "2") {
        
        my $cmd = "cat";
        foreach my $core_name (@core_names) {
            
            my $file = "$core_name" . "_$end.fastq";
            $cmd .= " $file ";
            push (@tmp_files, $file);
        }
        $cmd .= ">$prefix" . "_$end.fastq";
        &process_cmd($cmd);
    }
    
    ## remove tmp files:
    foreach my $file (@tmp_files) {
        unlink($file);
    }
}


exit(0);


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
