package COMMON;

use strict;
use warnings;
use Carp;

$ENV{LC_ALL} = 'C'; # needed for sorting order.

####
sub get_sort_exec {
    my ($num_threads) = @_;

    # check it like so:
    #  perl -MCOMMON -e 'print COMMON::get_sort_exec(4);'

    my $sort_exec = `which sort`;
    unless ($sort_exec =~ /\w/) {
        confess "Error, cannot find sort utility";
    }
    $sort_exec =~ s/\s//g;
    
    my $help_text = `$sort_exec --help`;  
    if ($help_text =~ m|--parallel|) {
        ## could do simple versioning check, but I don't remember which version started using parallel
        $sort_exec = "$sort_exec --parallel=$num_threads";
    }
    
    return($sort_exec);
}

1; #EOM
