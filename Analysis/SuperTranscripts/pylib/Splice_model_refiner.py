#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

import TGraph
import TNode
import Node_alignment

from Compact_graph_whole import Compact_graph_whole
from Compact_graph_partial import Compact_graph_partial

logger = logging.getLogger(__name__)


MAX_MM_RATE = 0.05 


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

    graph_compactor = Compact_graph_whole()
    graph_compactor.compact_unbranched(refined_tgraph)

    refined_tgraph.draw_graph("ladeda.linear_compact.dot")
    
    ###########
    #
    
    for allowed_variants in (0, 1, 2):
        graph_compactor.compact_graph(refined_tgraph, allowed_variants)
        refined_tgraph.draw_graph("ladeda.compact.m{}.dot".format(allowed_variants))
        

    ## now extract prefix and suffix matches
    #compact_graph_partial(refined_tgraph)
    

    refined_tgraph.draw_graph("ladeda.final.dot")
    
    # convert compacted graph into a node alignment obj
    return(node_alignment_obj)


