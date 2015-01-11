#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 rsem.file eXpress.file count|FPKM\n\n";

my $rsem_file = $ARGV[0] or die $usage;
my $express_file = $ARGV[1] or die $usage;
my $method = $ARGV[2] or die $usage;

unless ($method =~ /^(count|FPKM)$/) { 
    die $usage;
}

main: {

    my %data;
    &add_RSEM_data(\%data, $rsem_file, $method);
    &add_eXpress_data(\%data, $express_file, $method);


    print join("\t", "", "RSEM", "eXpress") . "\n";
    foreach my $id (keys %data) {
        my $rsem_val = $data{$id}->{RSEM};
        unless (defined $rsem_val) {
            $rsem_val = -1;
        }
        my $express_val = $data{$id}->{eXpress};
        unless (defined $express_val) {
            $express_val = -1;
        }

        if ($rsem_val > 0 && $express_val > 0) {
            $rsem_val = "NA" if $rsem_val < 0;
            $express_val = "NA" if $express_val < 0;
            
            print join("\t", $id, $rsem_val, $express_val) . "\n";
        }
    }
    
    exit(0);
}

####
sub add_RSEM_data {
    my ($data_href, $rsem_file, $method) = @_;

    my $data_field = ($method eq 'count') ? 4 : 6;
   
    open (my $fh, $rsem_file) or die "Error, cannot open file $rsem_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $id = $x[0];
        my $val = $x[$data_field];

        $data_href->{$id}->{RSEM} = $val;
    }
    close $fh;

    return;
}

####
sub add_eXpress_data {
    my ($data_href, $express_file, $method) = @_;

    my $data_field = ($method eq 'count') ? 7 : 10;

    open (my $fh, $express_file) or die "Error, cannot open file $express_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $id = $x[1];
        my $val = $x[$data_field];

        $data_href->{$id}->{eXpress} = $val;
    }
    close $fh;

    return;
}

        
