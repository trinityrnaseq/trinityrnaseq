#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 iworm_bundles_fasta_file sorted_reads_to_components_file\n\n\n";

my $iworm_bundles_fasta_file = $ARGV[0] or die $usage;
my $reads_to_components_file = $ARGV[1] or die $usage;


main: {

    my %component_to_LRs;
    {
        open (my $fh, $reads_to_components_file) or die "Error, cannot open file: $reads_to_components_file";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $component_id = $x[0];
            my $read_name = $x[1];

            if ($read_name =~ /^>LR\$\|/) {
                # got a long read!
                my $seq = $x[3];
                $component_to_LRs{$component_id} .= "X$seq";
            }
        }
        close $fh;
    }

    ## now add them to the compoennts:
    open (my $fh, $iworm_bundles_fasta_file) or die "Error, cannot open file: $iworm_bundles_fasta_file";
    while (my $header = <$fh>) {
        
        chomp $header;
        
        my $seq = <$fh>;

        chomp $seq;

        $header =~ /^>s_(\d+)/ or die "Error, cannot decipher component id";
        my $component_id = $1;
        
        if (my $LRs = $component_to_LRs{$component_id}) {
            $seq .= $LRs;
        }

        print "$header\n$seq\n";
    }

    exit(0);
}


