#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

from TGraph import *
from TNode import *
from Node_alignment import *


logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def refine_alignment(node_alignment_obj):

    """
    Create a new splice graph based on the node alignment obj.

    Since some nodes may show up as duplicate (repeat) nodes, assign each a unique ID
    """

    transcript_names = node_alignment_obj.get_transcript_names()
    aligned_nodes = node_alignment_obj.get_aligned_nodes()

    width = node_alignment_obj.width()

    new_node_list = list()
    orig_node_list = list()

    refined_tgraph = TGraph("^^SGRAPH2^^")
    
    for i in range(0,width):
        repr_node = node_alignment_obj.get_representative_column_node(i)
        orig_node_list.append(repr_node)

        transcripts = repr_node.get_transcripts()
        
        new_node = refined_tgraph.get_node(transcripts, "loc_" + str(i), repr_node.get_seq())
        new_node_list.append(new_node)

    # build graph
    for iso_node_alignment in aligned_nodes:

        prev = None
        for i in range(0,width):

            if iso_node_alignment[i] != None:
                if prev != None:
                    prev.add_next_node(new_node_list[i])
                prev = new_node_list[i]
            
            

    for node in new_node_list:
        print(node.toString() + " " + node.get_seq())


    refined_tgraph.draw_graph("ladeda.dot")
    


    return(node_alignment_obj)

    
    
