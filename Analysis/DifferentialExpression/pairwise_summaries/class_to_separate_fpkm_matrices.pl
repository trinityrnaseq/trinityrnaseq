#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 matrix.fpkm venn.class\n\n";

my $fpkm_matrix = $ARGV[0] or die $usage;
my $venn_class = $ARGV[1] or die $usage;

  
