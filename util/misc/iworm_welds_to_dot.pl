#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 iworm_cluster_welds_graph.txt.sorted.wIwormNames\n\n";

my $iworm_welds = $ARGV[0] or die $usage;

main: {

    open (my $fh, $iworm_welds) or die "Error, cannot open file $iworm_welds";
    
    print "digraph G {\n";

    while (<$fh>) {
        chomp;
        my @x = split(/\s+/);
        my $iworm_A = $x[1];
        my $iworm_B = $x[4];
        $iworm_A =~ s/;/_/;
        $iworm_B =~ s/;/_/;
        print "    $iworm_A->$iworm_B\n";
        
    }
    print "}\n";
    
    exit(0);
}


        
