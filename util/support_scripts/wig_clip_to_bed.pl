#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use WigParser;
use Gene_obj;


my $usage = "usage: $0 clip.wig\n\n";

my $clip_wig = $ARGV[0] or die $usage;

main: {

    
    my %mol_to_clips = &get_clip_pts($clip_wig);

    foreach my $mol (keys %mol_to_clips) {
        
        my @clip_pos = @{$mol_to_clips{$mol}};
        
        foreach my $pos (@clip_pos) {
            
            my ($end5, $end3) = ($pos-1, $pos+1);
            
            my $gene_obj = new Gene_obj();
            $gene_obj->populate_gene_object( {$end5=>$end3}, {$end5=>$end3});

            $gene_obj->{asmbl_id} = $mol;
            
            $gene_obj->{com_name} = "$mol-C:$pos";
            
            print $gene_obj->to_BED_format();
        }
    }

    exit(0);
}


####
sub get_clip_pts {
    my ($jaccard_clips) = @_;

    my %trans_to_clips;

    my $trans_acc = "";

    open (my $fh, $jaccard_clips) or die "Error, cannot open file $jaccard_clips";
    while (<$fh>) {
        chomp;
        if (/^variableStep chrom=(\S+)/) {
            $trans_acc = $1;
        }
        elsif (/^(\d+)/) {
            my ($coord, $val) = split(/\t/);
            if ($val) {
                push (@{$trans_to_clips{$trans_acc}}, $coord);
            }
        }
    }
    close $fh;

    return(%trans_to_clips);
}


        
