#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

import TGraph
import TNode
import Node_alignment
import Topological_sort

from Compact_graph_whole import Compact_graph_whole
from Compact_graph_partial import Compact_graph_partial
import TGLOBALS

logger = logging.getLogger(__name__)


MAX_MM_RATE = 0.05 


def refine_alignment(node_alignment_obj, reset_node_ids=False):

    """
    Create a new splice graph based on the node alignment obj.

    Since some nodes may show up as duplicate (repeat) nodes, assign each a unique ID
    """
    
    logger.debug("refine_alignment({})".format(node_alignment_obj))


    # convert to splice graph
    refined_tgraph = node_alignment_obj.to_splice_graph("^^SGRAPH2^^", reset_node_ids)

    if TGLOBALS.DEBUG:
        refined_tgraph.draw_graph("ladeda.pre.dot")
    
    graph_compactor = Compact_graph_whole()
    graph_compactor.compact_unbranched(refined_tgraph)

    if TGLOBALS.DEBUG:
        refined_tgraph.draw_graph("ladeda.linear_compact.dot")
    
    ###########
    #
    
    for allowed_variants in (0, 1, 2):
        graph_compactor.compact_graph(refined_tgraph, allowed_variants)
        if TGLOBALS.DEBUG:
            refined_tgraph.draw_graph("ladeda.compact.m{}.dot".format(allowed_variants))
        

    ## now extract prefix and suffix matches

    partial_graph_compactor = Compact_graph_partial()
    partial_graph_compactor.compact_graph(refined_tgraph)
    
    if TGLOBALS.DEBUG:
        refined_tgraph.draw_graph("ladeda.final.dot")
    
    # convert compacted graph into a node alignment obj

    splice_graph_node_alignment = splice_graph_to_node_alignment(refined_tgraph)

    splice_graph_node_alignment = remove_redundant_paths(splice_graph_node_alignment)

    
    return(splice_graph_node_alignment)



def splice_graph_to_node_alignment(tgraph):

    topologically_sorted_nodes = Topological_sort.Topological_sort.topologically_sort(tgraph.get_all_nodes())

    logger.debug("Topologically sorted nodes: " + str(topologically_sorted_nodes))

    # index loc node ids
    aligned_loc_id_pos = dict()
    for i in range(0, len(topologically_sorted_nodes)):
        loc_id = topologically_sorted_nodes[i].get_loc_id()
        aligned_loc_id_pos[loc_id] = i


    new_alignments = list()
    transcript_ids = set()
    for node in topologically_sorted_nodes:
        transcript_ids = transcript_ids.union(node.get_transcripts())

    transcript_ids = list(transcript_ids)
    
    for transcript_id in transcript_ids:
        new_alignment = [None for i in topologically_sorted_nodes]
        for node in topologically_sorted_nodes:
            if transcript_id in node.get_transcripts():
                loc_id = node.get_loc_id()
                new_idx = aligned_loc_id_pos[loc_id]
                new_alignment[new_idx] = node
        new_alignments.append(new_alignment)

    splice_graph_node_alignment = Node_alignment.Node_alignment(tgraph.get_gene_id(), transcript_ids, new_alignments)

    logger.debug("Splice graph node_alignment: " + str(splice_graph_node_alignment))

    return(splice_graph_node_alignment)



def remove_redundant_paths(node_alignment):

    transcript_names = node_alignment.get_transcript_names()
    aligned_nodes = node_alignment.get_aligned_nodes()

    num_transcripts_before_reduction = len(transcript_names)

    # do all pairwise comparisons
    # check for containments

    containments = set()

    for i in range(0,len(aligned_nodes)-1):
        for j in range(i+1, len(aligned_nodes)):
            if a_contains_b(aligned_nodes[i], aligned_nodes[j]):
                containments.add(j)
            elif a_contains_b(aligned_nodes[j], aligned_nodes[i]):
                containments.add(i)

    if containments:
        adj_transcript_names = list()
        adj_aligned_nodes = list()
        for i in range(0, len(aligned_nodes)):
            if i not in containments:
                adj_transcript_names.append(transcript_names[i])
                adj_aligned_nodes.append(aligned_nodes[i])

        adj_splice_graph_node_alignment = Node_alignment.Node_alignment(adj_transcript_names, adj_aligned_nodes)

        num_after_reduction = len(adj_transcript_names)
        
        logger.debug("Containments found, reporting reduced set {} of {} = {:.2f}%".format(
            num_after_reduction, num_transcripts_before_reduction,
            num_after_reduction/num_transcripts_before_reduction*100))
        
        return adj_splice_graph_node_alignment

    else:
        logger.debug("No containments found")
        return node_alignment # no changes


def get_first_node_idx(node_list):

    # find starting place for comparison
    begin_idx = -1
    for i in range(0,len(node_list)):
        if node_list[i] is not None:
            begin_idx = i
            break

    if begin_idx < 0:
        raise RuntimeError("Error, didn't find first non-none value among {} and {}".format(node_list_A, node_list_B))

    return begin_idx


def get_end_node_idx(node_list):
    
    # find ending place for comparison
    end_idx = -1

    for i in reversed(range(0, len(node_list))):
        if node_list[i] is not None:
            end_idx = i
            break
    
    if end_idx < 0:
        raise RuntimeError("Error, didn't find last non-none value among {} and {}".format(node_list_A, node_list_B))

    return end_idx



def a_contains_b(node_list_A, node_list_B):

    A_start = get_first_node_idx(node_list_A)
    A_end = get_end_node_idx(node_list_A)

    B_start = get_first_node_idx(node_list_B)
    B_end = get_end_node_idx(node_list_B)

    if not (A_start <= B_start and A_end >= B_end):
        return False # no containment

    # ensure that in overlapping region, they have identical nodes
    for i in range(B_start, B_end+1):
        # if we see any difference, then not compatible 
        if node_list_A[i] != node_list_B[i]:
            return False

    # must be compatible and contained
    return True
