#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use lib ("$FindBin::Bin/../../PerlLib", "$FindBin::Bin/PerlLib");
use Gene_obj;
use GFF3_utils;

use SegmentGraph;

my $usage = "usage: $0 genes.gff3\n\n";

my $gff3_file = $ARGV[0] or die $usage;

main: {


    my $gene_obj_indexer_href = {};
    ## associate gene identifiers with contig id's.
    my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

    foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
        
        my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
        foreach my $gene_id (@gene_ids) {
            my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
                
            my @additional_isoforms = $gene_obj_ref->get_additional_isoforms();
            
            unless (@additional_isoforms) {
                next;
            }

            my $segment_graph = new SegmentGraph();
                    
            foreach my $isoform ($gene_obj_ref, @additional_isoforms) {

                my $isoform_id = $isoform->{Model_feat_name};
                
            
                my @exons = $isoform->get_exons();
                
                foreach my $exon (@exons) {
                    
                    my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
                    
                    $segment_graph->add_segment($lend, $rend, $isoform_id);

                }
            }
            
            
            #print "Graph for: $gene_id\n";
            #print $segment_graph->toString() . "\n\n\n";
            
            &identify_longest_path_between_two_sites_of_variation($gene_id, $segment_graph);
            
            

        }
        
        
 
    }

    exit(0);
}
                    
####
sub identify_longest_path_between_two_sites_of_variation {
    my ($gene_id, $segment_graph) = @_;

    my @isoforms = $segment_graph->identify_all_owners();

    print "$gene_id\tISOFORMS: @isoforms\n";
    
    for (my $i = 0; $i < $#isoforms; $i++) {

        my $isoform_i = $isoforms[$i];

        for (my $j = $i + 1; $j <= $#isoforms; $j++) {
            
            my $isoform_j = $isoforms[$j];
            
            my @longest_path = &find_longest_path_between_iso_pair($segment_graph, $isoform_i, $isoform_j);
            
            
            


        }

    }


    return;
}
    


####
sub find_longest_path_between_iso_pair {
    my ($segment_graph, $iso_A, $iso_B) = @_;
    
    my @nodes = $segment_graph->get_all_nodes();  ## ordered left to right.

    my %seen;
    
    my @long_path;
    my @all_long_paths;;
    
    
  path_search:
    foreach my $seed_node (@nodes) {
        my $node_ID = $seed_node->get_ID();

        if ($seen{$node_ID}) { next; }
        $seen{$node_ID} = 1;

        unless ($seed_node->has_owners($iso_A, $iso_B)) {
            next;
        }

        ## ensure left-branched with respect to A,B
        
        push (@long_path, $seed_node);
        my $next_node = $seed_node;
        my $forward_OK = 0;
        while (1) {
            my @next_nodes = $next_node->get_next_nodes();
            unless (@next_nodes) {
                last; # reached end w/o finding terminting node.
            }
          
            my $found_next_node = 0;
            my $found_A = 0;
            my $found_B = 0;
     
          next_node_candidate_search:
            foreach my $next_node_candidate (@next_nodes) {
                if ($next_node_candidate->has_owners($iso_A, $iso_B)) {
                    ## should only be one such path. No cycles allowed in gene structures.
                    
                    $found_next_node = 1;
                    push (@long_path, $next_node);
                    $next_node = $next_node_candidate;
                    last next_node_candidate_search;
                }
                elsif ($next_node_candidate->has_owners($iso_A)) {
                    $found_A = 1;
                }
                elsif ($next_node_candidate->has_owners($iso_B)) {
                    $found_B = 1;
                }
            }
            if (! $found_next_node) {
                if ($found_A && $found_B) {
                    ## excellent, found branched diff.
                    push (@all_long_paths, [@long_path]);
                    @long_path = ();
                                        
                }
                else {
                    ## doesn't meet our requirements for branched long path diff.
                    @long_path = (); # clear it out, no save.
                    
                }
                last; # break while.
            }
        }

    }
    
    
}
