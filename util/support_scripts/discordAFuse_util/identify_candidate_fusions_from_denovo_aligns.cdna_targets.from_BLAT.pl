#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;

use FindBin;
use lib "$FindBin::Bin/../../../PerlLib";
require "overlapping_nucs.ph";

my $MAX_ALLOWED_PCT_OVERLAP = 20;


my $usage = "\n\n\tusage: $0 blat.gff3 cdna_targets.fasta\n\n";

my $blat_gff3 = $ARGV[0] or die $usage;
my $cdna_targets_fasta = $ARGV[1] or die $usage;

main: {


    my $headers_file = "$cdna_targets_fasta.headers";
    unless (-s $headers_file) {
        die "Error, cannot locate $headers_file";
    }

    my %trans_id_to_gene_id = &get_trans_to_gene_mapping($headers_file);
     
    my %trans_to_gene_alignments = &get_trans_to_gene_alignments($blat_gff3, \%trans_id_to_gene_id);
    
    my %tiered_n_merged_gene_alignments = &tier_and_merge_best_alignments(\%trans_to_gene_alignments);
    
    foreach my $trans_id (keys %tiered_n_merged_gene_alignments) {
        
        my @hits = @{$tiered_n_merged_gene_alignments{$trans_id}};
        @hits = sort {$a->{trans_lend}<=>$b->{trans_lend}} @hits;

        if (scalar @hits < 2) {
            # shouldn't happen as filtering is done before this.
            die "Error, have fewer than 2 candidate fusion parts";
        }
        
        ## Report fusion pairs:
        my $geneA = shift @hits;
        while (@hits) {
            my $geneB = shift @hits;

            if ($geneA->{gene_orient} eq $geneB->{gene_orient}) {

                my ($left_gene, $right_gene) = ($geneA->{gene_orient} eq '+') ? ($geneA, $geneB) : ($geneB, $geneA);
                
                my $fusion_name = join("--", $left_gene->{gene_name}, $right_gene->{gene_name});

                print join("\t", $fusion_name, $trans_id,
                           $left_gene->{gene_name}, $left_gene->{trans_lend}, $left_gene->{trans_rend}, $left_gene->{gene_orient}, $left_gene->{per_id} . "\%ID",
                           $right_gene->{gene_name}, $right_gene->{trans_lend}, $right_gene->{trans_rend}, $right_gene->{gene_orient}, $right_gene->{per_id} . "\%ID") . "\n";

            }
            
            $geneA = $geneB;
        }
    }
    
    exit(0);


}

####
sub tier_and_merge_best_alignments {
    my ($trans_to_gene_alignments_href) = @_;

    my %trans_id_to_tiers;

    foreach my $trans_id (keys %$trans_to_gene_alignments_href) {
        
        my @hits = @{$trans_to_gene_alignments_href->{$trans_id}};
        
        my @merged_gene_hits = &merge_hits_by_gene(@hits);

        unless (scalar @merged_gene_hits > 1) {
            next;
        }

        foreach my $hit (@merged_gene_hits) {
            $hit->{score} = ($hit->{trans_rend} - $hit->{trans_lend} + 1) * $hit->{per_id};
        }

        @hits = reverse sort {$a->{score}<=>$b->{score}} @merged_gene_hits;
        
        
        my @tiered_hits;
        foreach my $hit (@hits) {
            if (&insufficient_overlap_existing_tier(\@tiered_hits, $hit)) {
                push (@tiered_hits, $hit);
            }
        }
        if (scalar @tiered_hits > 1) {
            $trans_id_to_tiers{$trans_id} = [@tiered_hits];
        }
    }

    return(%trans_id_to_tiers);
    
}


####
sub insufficient_overlap_existing_tier {
    my ($tiered_aref, $hit) = @_;

    foreach my $tier_entry (@$tiered_aref) {
        
        my $tier_gene_name = $tier_entry->{gene_name};
        my $tier_lend = $tier_entry->{trans_lend};
        my $tier_rend = $tier_entry->{trans_rend};
        my $tier_seg_len = $tier_rend - $tier_lend + 1;


        my $hit_gene_name = $hit->{gene_name};
        my $hit_lend = $hit->{trans_lend};
        my $hit_rend = $hit->{trans_rend};
        my $hit_len = $hit_rend - $hit_lend + 1;


        #print "TESTING $tier_gene_name vs. $hit_gene_name : ($tier_lend-$tier_rend) vs. ($hit_lend, $hit_rend)\n";

        
        if (&coordsets_overlap([$tier_lend, $tier_rend], [$hit_lend, $hit_rend])) {

            my $smaller_len = ($tier_seg_len < $hit_len) ? $tier_seg_len : $hit_len;

            my $pct_overlap = &nucs_in_common($tier_lend, $tier_rend, $hit_lend, $hit_rend) / $smaller_len * 100;
            
            #print "($tier_lend-$tier_rend) vs. ($hit_lend, $hit_rend)\tPCT_Overlap: $pct_overlap\n";

            if ($pct_overlap > $MAX_ALLOWED_PCT_OVERLAP) {
                return(0);
            }
        }
        else {
            #print "($tier_lend-$tier_rend) vs. ($hit_lend, $hit_rend)\tNO_overlap\n";
        }

    }
    return(1); # no sufficient overlap found
    
}
        
        



####
sub merge_hits_by_gene {
    my @hits = @_;

    my %gene_to_hit_list;

    foreach my $hit (@hits) {
        my $gene_name = $hit->{gene_name};
        
        push (@{$gene_to_hit_list{$gene_name}}, $hit);
    }


    ## merge them
    my @merged_gene_hits;

    foreach my $gene_name (keys %gene_to_hit_list) {
        my @hits = @{$gene_to_hit_list{$gene_name}};
        
        ## todo: check for actual overlaps
        ## first, keeping it super easy and selecting range of matches 

        my @lend_coords;
        my @rend_coords;
        
        my $sum_len = 0;
        my $sum_per_id_n_len = 0;

        foreach my $hit (@hits) {
            push (@lend_coords, $hit->{trans_lend});
            push (@rend_coords, $hit->{trans_rend});
        
            my $seg_len = $hit->{trans_rend} - $hit->{trans_lend} + 1;
            my $per_id = $hit->{per_id};

            $sum_per_id_n_len += $seg_len * $per_id;

            $sum_len += $seg_len;
        }

        @lend_coords = sort {$a<=>$b} @lend_coords;
        @rend_coords = sort {$a<=>$b} @rend_coords;

        my $range_lend = shift @lend_coords;
        my $range_rend = pop @rend_coords;


        my $avg_per_id = sprintf("%.1f", $sum_per_id_n_len / $sum_len);
        my $hit = shift @hits;
        $hit->{trans_lend} = $range_lend;
        $hit->{trans_rend} = $range_rend;
        $hit->{per_id} = $avg_per_id;

        push (@merged_gene_hits, $hit);
    }

    return(@merged_gene_hits);
}



####
sub get_trans_to_gene_alignments {
    my ($blat_gff3, $trans_id_to_gene_id_href) = @_;

    my %trans_to_gene_hits;

    open (my $fh, $blat_gff3) or die "Error, cannot open file $blat_gff3";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $target_trans_id = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $per_id = $x[5];
        my $orient = $x[6];
        my $info = $x[8];

        $info =~ /Target=(\S+) (\d+) (\d+)/ or die "Error, cannot parse hit info from $info";
        
        my $trans_id = $1;
        my $trans_end5 = $2;
        my $trans_end3 = $3;
        

        my $gene_name = $trans_id_to_gene_id_href->{$target_trans_id} or die "Error, no gene name for $target_trans_id";
        
        push (@{$trans_to_gene_hits{$trans_id}}, { gene_name => $gene_name,
                                                   trans_lend => $trans_end5,
                                                   trans_rend => $trans_end3,
                                                   gene_orient => $orient,
                                                   per_id => $per_id,
              }
            );

    }

    return(%trans_to_gene_hits);
}


####
sub get_trans_to_gene_mapping {
    my ($headers_file) = @_;
    
    my %trans_id_to_gene_id;

    open (my $fh, $headers_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($trans_id, $gene_id, $gene_name) = split(/\s+/);
        
        if ($gene_name) {
            $trans_id_to_gene_id{$trans_id} = $gene_name;
        }
        else {
            $trans_id_to_gene_id{$trans_id} = $gene_id;
        }

    }
    close $fh;


    return(%trans_id_to_gene_id);
}
