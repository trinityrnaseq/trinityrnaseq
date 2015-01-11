#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 --labels sampleA sampleB [...] --RSEM_counts sampleA.RSEM.isoform.results sampleB.RSEM.isoform.results [...]\n\n";

unless (@ARGV) {
    die $usage;
}

unless (scalar @ARGV >= 6) {
    die $usage;
}


my @sample_labels;
my @rsem_files;

my $list_aref;
while (@ARGV) {
    my $param = shift;
    if ($param eq "--labels") {
        $list_aref = \@sample_labels;
    }
    elsif ($param eq "--RSEM_counts") {
        $list_aref = \@rsem_files;
    }
    else {
        unless ($list_aref) {
            die "Error, cannot determine input category for entry $param";
        }
        push (@$list_aref, $param);
    }
}

unless (scalar(@sample_labels) == scalar(@rsem_files)) {
    die "Error, the number of sample labels != number of rsem files:\n"
        . "sample labels: " . join(", ", @sample_labels) . "\n"
        . "count files: " . join(", ", @rsem_files) . "\n";
}

unless (scalar(@sample_labels) >= 2 && scalar(@rsem_files) >= 2) {
    die "numbers of sample labels and rsem files do not match";
}

main: {


    my %data;
    
    foreach my $file (@rsem_files) {
        
        open (my $fh, $file) or die "Error, cannot open file $file";
        while (<$fh>) {
            chomp;
            my ($acc, $count, $rest) = split(/\t/);
            $data{$acc}->{$file} = $count;
        }
        close $fh;
    }
    
    print join("\t", "", @sample_labels) . "\n";
    foreach my $acc (keys %data) {
        
        print "$acc";
        
        foreach my $file (@rsem_files) {
            
            my $count = $data{$acc}->{$file};
            unless (defined $count) {
                $count = "NA";
            }
            
            print "\t$count";
            
        }
        
        print "\n";
        
    }
    
    
    exit(0);
}
    
        
