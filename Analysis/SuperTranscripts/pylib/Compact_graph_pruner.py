#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

import TGraph
import TNode
import TGLOBALS


logger = logging.getLogger(__name__)



class Compact_graph_pruner:

    def __init__(self):

        pass


    def remove_burrs(self, tgraph, max_burr_length):
        nodes = tgraph.get_all_nodes()

                
        for node in nodes:
            if node.is_burr() and len(node.get_seq()) <= max_burr_length:
                logger.debug("Compact_graph_pruner.remove_burrs - pruning burr node: {}".format(node))
                tgraph.prune_node(node)
    
    
    
    def pop_small_bubbles(self, tgraph, max_bubble_node_length):
        """
        bubble structure:

             X1
            /  \
        ---A    B--
            \  /
             X2


        keep X1 or X2, ideally based on which is in a more highly expressed isoform structure.

        A node in a bubble is defined as:
        - has single parent and single child
        - the parents of child(A) and the children of the parent(B) are identical. (simple bubble structure)
        
        """

        bubble_node_lists = self._get_bubbles(tgraph, max_bubble_node_length)

        for bubble_node_list in bubble_node_lists:
            # TODO: should sort by expression, select representative node based on highest expr
            repr_node = bubble_node_list.pop()
            for sister_bubble_node in bubble_node_list:
                repr_node.add_transcripts(sister_bubble_node.get_transcripts())
                tgraph.prune_node(sister_bubble_node)





    def _get_bubbles(self, tgraph, max_bubble_node_length):

        bubble_node_lists = list()

        node_found_in_bubble = dict()  # store T if in bubble

        for node in tgraph.get_all_nodes():
            
            if node in node_found_in_bubble:
                # already identified as within a bubble as seeded by a different node.
                continue

            if len(node.get_prev_nodes()) != 1:
                continue
            if len(node.get_next_nodes()) != 1:
                continue

            child_node_A = node.get_prev_nodes().pop()
            parent_node_B = node.get_next_nodes().pop()

            bubble_nodes = child_node_A.get_next_nodes() & parent_node_B.get_prev_nodes()
            if len(bubble_nodes) < 2:
                # not a bubble
                continue

            # potential bubble...  check to see that each has only one parent and child
            selected_bubble_nodes = list()
            for bubble_node in bubble_nodes:
                if len(bubble_node.get_prev_nodes()) != 1:
                    continue
                if len(bubble_node.get_next_nodes()) != 1:
                    continue
                if len(bubble_node.get_seq()) > max_bubble_node_length:
                    continue
                
                selected_bubble_nodes.append(bubble_node)

            if len(selected_bubble_nodes) < 2:
                # no simple bubble here
                continue

            # got a simple bubble.  store it.
            bubble_node_lists.append(selected_bubble_nodes)
            logger.debug("Found bubble: {} -- {} -- {}".format(child_node_A, selected_bubble_nodes, parent_node_B))


            for selected_bubble_node in selected_bubble_nodes:
                node_found_in_bubble[selected_bubble_node] = True

        return bubble_node_lists


                         
        

        
        
