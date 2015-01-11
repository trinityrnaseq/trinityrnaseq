#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 coordSorted.sam chyrsalis_dir/\n\n";

my $coordSorted_sam_file = $ARGV[0] or die $usage;
my $chrysalis_dir = $ARGV[1] or die $usage;

my $iworm_graph_file = "$chrysalis_dir/GraphFromIwormFasta.out";
my $component_listing_file = "iworm_bundle_file_listing.txt";

main: {

    ## get link between components and iworm assemblies
    my %iworm_acc_to_component = &get_iworm_to_component_mappings($iworm_graph_file);

    my %component_file_prefix = &parse_component_listing($component_listing_file);
    
    my %seen_component;
    my $prev_component = "";
    my $ofh;


    my $counter = 0;
    
    ## write the reads
    open (my $sam_fh, $coordSorted_sam_file) or die $!;
    while (<$sam_fh>) {
        $counter++;

        if ($counter % 1000 == 0) {
            print STDERR "\r[$counter]     ";
        }
        

        my @x = split(/\t/);
        my ($read_acc, $iworm_acc, $sequence) = ($x[0], $x[2], $x[9]);

        my $component_id = $iworm_acc_to_component{$iworm_acc};
        unless (defined $component_id) {
            print STDERR "Error, iworm_acc: $iworm_acc was not assigned to a component.\n";
            next;
        }
        
        my $component_prefix = $component_file_prefix{$component_id};

        unless (defined $component_prefix) {
            #print STDERR "Error, no file prefix for component: $component_id\n";
            next;
        }
        
        if ($component_id ne $prev_component) {
            close $ofh if $ofh;
            
            my $reads_file = "$component_prefix.raw.fasta";
            
            if ($seen_component{$component_id}) {
                open ($ofh, ">>$reads_file") or die $!;
            }
            else {
                open ($ofh, ">$reads_file") or die $!;
                $seen_component{$component_id} = 1;
            }
        }
        if ($read_acc =~ /(-\d)$/) {
            $read_acc =~ s/-(\d)$/\/$1/;
        }
        
        print $ofh ">$read_acc\n$sequence\n";
        #print ">$read_acc\n$sequence\n";

        $prev_component = $component_id;
        
        
        
        
    }
    close $sam_fh;
    
    close $ofh;
    
    
    {
        open (my $ofh, ">$chrysalis_dir/readcounts.out") or die $!;
        print $ofh $counter . "\n";
        close $ofh;
    }

    
}

####
sub parse_component_listing {
    my ($component_listing_file) = @_;

    my %component_to_prefix;

    open (my $fh, $component_listing_file) or die $!;
    while (<$fh>) {
        chomp;
        my $path = $_;
        
        $path =~ /^(.*comp(\d+))/ or die "Error, cannot parse component info from $path";
        my $comp_prefix = $1;
        my $component_no = $2;
        $component_to_prefix{$component_no} = $1;
    }
    close $fh;

    return (%component_to_prefix);
}



####
sub get_iworm_to_component_mappings {
    my ($iworm_graph_file) = @_;

    my %iworm_to_comp;

    open (my $fh, "$iworm_graph_file") or die $!;
    while (<$fh>) {
        chomp;
        if (/^>Component_(\d+) .* \[iworm>(a\d+;\d+)/) {
            my $component_no = $1;
            my $iworm_acc = $2;
            $iworm_to_comp{$iworm_acc} = $component_no;
        }
    }
    close $fh;

    return(%iworm_to_comp);
}
