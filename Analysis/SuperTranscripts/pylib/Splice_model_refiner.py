#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

import TGraph
import TNode
import Node_alignment



logger = logging.getLogger(__name__)
#logger.addHandler(logging.NullHandler())

def refine_alignment(node_alignment_obj):

    """
    Create a new splice graph based on the node alignment obj.

    Since some nodes may show up as duplicate (repeat) nodes, assign each a unique ID
    """
    
    logger.debug("refine_alignment({})".format(node_alignment_obj))

    transcript_names = node_alignment_obj.get_transcript_names()
    aligned_nodes = node_alignment_obj.get_aligned_nodes()

    width = node_alignment_obj.width()

    new_node_list = list()
    orig_node_list = list()

    refined_tgraph = TGraph.TGraph("^^SGRAPH2^^")
    
    for i in range(0,width):
        repr_node = node_alignment_obj.get_representative_column_node(i)
        orig_node_list.append(repr_node)

        transcripts = repr_node.get_transcripts()
        
        new_node = refined_tgraph.get_node(transcripts, "loc_" + str(i), repr_node.get_seq())
        new_node_list.append(new_node)

    #############
    # build graph
    
    for iso_node_alignment in aligned_nodes:

        prev = None
        for i in range(0,width):

            if iso_node_alignment[i] != None:
                if prev != None:
                    refined_tgraph.add_edges([prev], [new_node_list[i]])
                prev = new_node_list[i]
            
            
    logger.debug("new graph node listing:")
    for node in new_node_list:
        logger.debug(node.toString() + " " + node.get_seq())


    refined_tgraph.draw_graph("ladeda.pre.dot")


    ###########
    #

    for allowed_variants in (0, 1):
        compact_graph(refined_tgraph, allowed_variants)
        refined_tgraph.draw_graph("ladeda.compact.m{}.dot".format(allowed_variants))
        
    # convert compacted graph into a node alignment obj


    return(node_alignment_obj)

    
def compact_graph(tgraph, num_allowed_variants):

    compacted_flag = True

    while compacted_flag:

        tgraph.clear_touch_settings() # start fresh
        compacted_flag = False # reset for this round

        for node in tgraph.get_all_nodes():

            if node.get_touched_val() > 0:
                continue

            # try compact upward
            prev_nodes = node.get_prev_nodes()
            if len(prev_nodes) > 1 and untouched(prev_nodes):
                if compact_upward(prev_nodes, num_allowed_variants):
                    compacted_flag = True

            next_nodes = node.get_next_nodes()
            if len(next_nodes) > 1 and untouched(next_nodes):
                if compact_downward(next_nodes, num_allowed_variants):
                    compacted_flag = True


def untouched(node_list):
    for node in node_list:
        if node.get_touched_val() > 1:
            return False

    return True


####################
## Upward Compaction

def compact_upward(prev_nodes, num_allowed_variants):

    logger.debug("compact_upward: {} max_allowed_variants: {}".format(prev_nodes, num_allowed_variants))

    prev_nodes = list(prev_nodes)

    for i in range(0, len(prev_nodes)-1):
        for j in range(i+1, len(prev_nodes)):
            if compact_upward_node_pair(prev_nodes[i], prev_nodes[j], num_allowed_variants):
                return True

    return False


def compact_upward_node_pair(node_A, node_B, num_allowed_variants):

    # switch A,B so A is the shorter sequence node

    if len(node_B.get_seq()) < len(node_A.get_seq()):
        (node_A, node_B) = (node_B, node_A)
    
    seqA = node_A.get_seq()
    seqB = node_B.get_seq()

    shorter_len = len(seqA)


    # reverse the strings and check number of mismatches
    seqA_rev = seqA[::-1]
    seqB_rev = seqB[::-1]

    num_mm = 0
    for i in range(0, shorter_len):
        if seqA_rev[i] != seqB_rev[i]:
            num_mm += 1

    if num_mm > num_allowed_variants:
        logger.debug("num_mm: {} exceeds {} for node pair: {} and {}".format(num_mm, num_allowed_variants, node_A, node_B))
        return False

    # go ahead and merge the two in their region of common overlap
    tgraph = node_A.get_graph()

    ##### operations needed:
    #
    #         A
    #           \
    #            -- X
    #           / 
    #         B
    #
    #   -add B transcripts to A
    #   -adjust all B-downward edges to A
    #   -if B is longer than A, set downard to A and trim B part to unique region
    #      otherwise, remove B altogether.
    #

    node_A.add_transcripts(node_B.get_transcripts())
    tgraph.add_edges( [node_A], node_B.get_next_nodes() )
    tgraph.prune_edges( [node_B], node_B.get_next_nodes() )

    node_A.touch()
    node_B.touch()


    if len(seqB) > shorter_len:
        delta = len(seqB) - shorter_len
        seqB = seqB[0:delta]
        node_B.set_seq(seqB)
        tgraph.add_edges([node_B], [node_A])
    else:
        # remove node_B prev edges and attach them to A
        tgraph.add_edges(node_B.get_prev_nodes(), [node_A])
        # remove B altogether
        tgraph.prune_node(node_B)

    logger.debug("compacted nodes: {} and {}".format(node_A, node_B))
    
    return True

######################
## Downward compaction

def compact_downward(next_nodes, num_allowed_variants):

    next_nodes = list(next_nodes)
    for i in range(0, len(next_nodes)-1):
        for j in range(i+1, len(next_nodes)):
            if compact_downward_node_pair(next_nodes[i], next_nodes[j], num_allowed_variants):
                return True

    return False


def compact_downward_node_pair(node_A, node_B, num_allowed_variants):

    return False

            
