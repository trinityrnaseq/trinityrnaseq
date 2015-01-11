#!/usr/local/bin/perl

package main;
our $SEE;

package GFF_maker;

use strict;
use warnings;
use Carp;
use Data::Dumper;

####
sub get_GFF_line {
    my $input_href = shift;

    ## field 1: Sequence Identifier
    my $seq_id = $input_href->{seq_id} or confess "need seq_id " . Dumper($input_href);

    ## field 2: Source (default '.')
    my $source = $input_href->{source} || '.';

    ## field 3: type (use SO:id syntax or corresponding token)
    ## -could do some better validation here.
    my $type = $input_href->{type};

    ## field 4: starting coordinate (lend)
    my $lend = $input_href->{lend};

    ## field 5: ending coordinate (rend)
    my $rend = $input_href->{rend};

    unless (defined $lend && defined $rend) {
        confess "lend or rend coordinate not defined: "  . Dumper($input_href);
    }

    if ($lend == 0) {
        print STDERR "Warning, lend coordinate is zero, changing it to base 1.\n";
        print STDERR Dumper($input_href);
        $lend = 1;
    }
    
    unless ($lend <= $rend) {
        confess " lend > rend, not allowed. " . Dumper($input_href);
    }
    
    ## field 6: score
    my $score = $input_href->{score} || '.';

    ## field 7: strand
    my $strand = $input_href->{strand};
    unless ($strand eq '+' || $strand eq '-') {
        confess "need strand " . Dumper ($input_href);
    }

    ## field 8: strand
    my $phase = $input_href->{phase};
    unless (defined $phase) {
        $phase = '.';
    }
    
    ## field 9: attributes
    my $attributes = $input_href->{attributes} or confess "need attributes " . Dumper($input_href);
    unless ($attributes =~ /ID=/) {
        confess "attributes requires ID field. " . Dumper($input_href);
    }

    my $textline = sprintf("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\n",
                           $seq_id,
                           $source,
                           $type,
                           $lend,
                           $rend,
                           $score,
                           $strand,
                           $phase,
                           $attributes);

    return ($textline);
    
}


    
    
    
1; #EOM

