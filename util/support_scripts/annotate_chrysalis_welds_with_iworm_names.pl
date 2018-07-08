#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 iworm.fasta iworm_cluster_welds_graph.txt\n\n";

my $iworm_fa = $ARGV[0] or die $usage;
my $iworm_cluster_welds_file = $ARGV[1] or die $usage;

 main: {


     my %counter_to_iworm;
     
     {
         my $counter = 0;
         open(my $fh, $iworm_fa) or die "Error, cannot open $iworm_fa file";
         while (<$fh>) {
             if (/^>(\S+)/) {
                 $counter_to_iworm{$counter} = $1;
                 $counter++;
             }
         }
         close $fh;
     }

     open(my $fh, $iworm_cluster_welds_file) or die "Error, cannot open file: $iworm_cluster_welds_file";
     while(<$fh>) {
         chomp;
         my ($idx_A, $ptr, $idx_B, @rest) = split(/\s+/);
         my $iworm_A = $counter_to_iworm{$idx_A};

         if (! defined $iworm_A) { die "Error, no iworm acc for index: $idx_A"; }
         
         my $iworm_B = $counter_to_iworm{$idx_B};

         if (! defined $iworm_B) { die "Error, no iworm acc for index: $idx_B"; }
         
         print join(" ", $idx_A, $iworm_A, $ptr, $idx_B, $iworm_B, @rest) . "\n";
     }
     close $fh;


}

exit(0);
     
