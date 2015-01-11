#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use GFF3_alignment_utils;
use WigParser;
use Data::Dumper;

my $usage = "usage: $0 gmap.gff3 jaccard.wig [lowest_in_window_size=300]\n\n";

my $gmap_gff3 = $ARGV[0] or die $usage;
my $jaccard_wig = $ARGV[1] or die $usage;
my $window_search = $ARGV[2] || 300;


main: {

    print STDERR "-indexing wig: $jaccard_wig\n";
    my $wig_parser = new WigParser($jaccard_wig);
    
    my $alignment_indexer_href = {};
    

    print STDERR "-parsing GFF3 file: $gmap_gff3\n";
    my %contig_to_alignment_ids = &GFF3_alignment_utils::index_GFF3_alignment_objs($gmap_gff3, $alignment_indexer_href);

    ## group multi-paths
    my %core_acc_to_alignments;


    print STDERR "-grouping alignments by acc.\n\n";
    foreach my $contig (keys %contig_to_alignment_ids) {

        my @align_ids = @{$contig_to_alignment_ids{$contig}};

        foreach my $align_id (@align_ids) {

            my $core_acc = $align_id;
            $core_acc =~ s/\.path\d+$//;

            push (@{$core_acc_to_alignments{$core_acc}}, $align_id);

        }
    }

    print STDERR "-examining chimeras.\n";
    foreach my $core_acc (keys %core_acc_to_alignments) {
        
        my @align_ids = @{$core_acc_to_alignments{$core_acc}};
        my $num_alignments = scalar(@align_ids);
        if ($num_alignments != 2) {
            next;
        }

        #print "Got: " . join(", ", @align_ids) . "\n";
        
        my @alignment_objs;
        my @mcoordsets;

        my @genome_map_entries;

        foreach my $align_id (@align_ids) {
            my $align_obj = $alignment_indexer_href->{$align_id} or die "Error, no alignment retrieved for $align_id";
            push (@alignment_objs, $align_obj);
        
            my @trans_coords = $align_obj->get_mcoords();
            my ($trans_left, $trans_right) = sort {$a<=>$b} @trans_coords;
            push (@mcoordsets, [$trans_left, $trans_right]);
        
            my $genome_acc = $align_obj->{genome_acc};
            
            my ($genome_lend, $genome_rend) = $align_obj->get_coords;
            push (@genome_map_entries, { 

                genome_acc => $genome_acc,
                genome_lend => $genome_lend,
                genome_rend => $genome_rend,
                trans_lend => $trans_coords[0],
                trans_rend => $trans_coords[1],
            } );
            
        }
        
        @mcoordsets = sort {$a->[0]<=>$b->[0]} @mcoordsets;

        #print "// $core_acc\n" . Dumper(\@mcoordsets) . "\n";
    

        my $clip_pt = $mcoordsets[0]->[1];

        my $wig_val = -1;
        my $num_single = -1;
        my $num_both = -1;
        eval {

            my @wig_array = $wig_parser->get_wig_array($core_acc, 1);
            
            if ($#wig_array > $clip_pt) {
                $wig_val = &find_lowest_in_window(\@wig_array, $clip_pt, $window_search); #$wig_array[$clip_pt];
            }
            if (! defined $wig_val) {
                $wig_val = -1;
            }
            if (ref $wig_val) {
                $num_single = $wig_val->[1];
                $num_both = $wig_val->[2];
                $wig_val = $wig_val->[0];
            }
            
        };
        if ($@) {
            print STDERR "No jaccard wig pt for $core_acc\n";
        }
        

       
        @genome_map_entries = sort {$a->{trans_lend}<=>$b->{trans_lend}} @genome_map_entries;
        
        my $out_text = join("\t", $core_acc, $clip_pt, $wig_val, $num_single, $num_both);
        foreach my $map_entry (@genome_map_entries) {
            $out_text .= "\t" . $map_entry->{genome_acc} . ":" 
                . $map_entry->{genome_lend} . "(" . $map_entry->{trans_lend} . ")"
                . "-"
                . $map_entry->{genome_rend} . "(" . $map_entry->{trans_rend} . ")";
        }

        print $out_text . "\n";
        
        
    }
        
    
    exit(0);






}


####
sub find_lowest_in_window {
    my ($wig_array_aref, $clip_pt, $window_search) = @_;

    my $left_search = int($clip_pt - ($window_search/2));
    if ($left_search < 1) { $left_search = 1; }

    my $right_search = int($clip_pt + ($window_search/2));
    if ($right_search > $#$wig_array_aref) {
        $right_search = $#$wig_array_aref;
    }

    my $wig_chosen;

    for (my $i = $left_search; $i <= $right_search; $i++) {
        
        my $wig_ref = $wig_array_aref->[$i];
        if (ref $wig_ref) {
            if ( (! defined $wig_chosen)
                 ||
                 $wig_ref->[0] < $wig_chosen->[0]
                 ||
                 ($wig_ref->[0] == $wig_chosen->[0] && $wig_ref->[2] > $wig_chosen->[2]) ) {

                $wig_chosen = $wig_ref;
            }
        }
    }

    return($wig_chosen);
}


    
