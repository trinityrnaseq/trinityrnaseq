#!/usr/bin/env perl

use strict;

while (<STDIN>) {
    tr/\t\n\000-\037\177-\377/\t\n  /d;
    print;
}
