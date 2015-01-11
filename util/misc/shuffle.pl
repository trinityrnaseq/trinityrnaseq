#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw (shuffle);

my @list;

while (<STDIN>) {
    push (@list, $_);
}

@list = shuffle @list;


foreach my $ele (@list) {
    print $ele;
}


exit(0);

