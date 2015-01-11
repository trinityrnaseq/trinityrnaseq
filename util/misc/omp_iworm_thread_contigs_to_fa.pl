#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 tmp.iworm.thread.txt\n\n";

my $file = $ARGV[0] or die $usage;

main: {
    open (my $fh, $file) or die $!;
    
    my $seq_counter = 0;
  reader:
    while (1) {
        
        my @vals;
        for (1..4) {
            my $line = <$fh>;
            if (eof($fh)) { last reader; }
            chomp $line;
            push (@vals, $line);
        }

        $seq_counter++;

        my $seq = pop @vals;
        print ">s$seq_counter\n$seq\n";

    }
    
    exit(0);
}

        
