#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging
import argparse
import collections
import numpy
import time

import TNode

logger = logging.getLogger(__name__)


class TGraph:

 
    def __init__(self, gene_id):

        self.node_cache = dict()
        self.gene_id = gene_id

    
    def get_node(self, transcript_id, loc_node_id, node_seq):
        
        """
        Instantiates Node objects, and stores them in a graph.

        *** use this method for instantiating all Node objects ***

        use clear_node_cache() to clear the graph
                
        """
                
        logger.debug("{}\t{}".format(loc_node_id, node_seq))
        
        if len(node_seq) == 0:
            raise RuntimeError("Error, non-zero length node_seq required for parameter")


        if loc_node_id in self.node_cache:
            node_obj = self.node_cache[ loc_node_id ]
            node_obj.add_transcripts(transcript_id)
            if node_obj.seq != node_seq:
                ## FIXME: when internal node is found as a terminal node, we're using the k-1mer now for the terminal node too.
                logger.debug("WARNING: have conflicting node sequences for {} node_id: {}\n".format(self.get_gene_id(),
                                                                                                    loc_node_id) +
                             "{}\n vs. \n{}".format(node_obj.seq, node_seq))
                if len(node_obj.seq) < len(node_seq) and re.search("{}$".format(node_obj.seq), node_seq):
                    return node_obj
                elif len(node_obj.seq) > len(node_seq) and re.search("{}$".format(node_seq), node_obj.seq):
                    node_obj.seq = node_seq # reset to shorter sequence, should be k-1 shorter
                    return node_obj
                else:
                    raise RuntimeError("Error, have conflicting node sequences for node_id: {}\n".format(node_id) +
                                   "{}\n vs. \n{}".format(node_obj.seq, node_seq))
            else:
                return node_obj
            
        else:
            # instantiate a new one
            node_obj = TNode.TNode(self, transcript_id, loc_node_id, node_seq)
            self.node_cache[ loc_node_id ] = node_obj
            return node_obj



    def get_all_nodes(self):
        return self.node_cache.values()
    
    def clear_node_cache(self):
        """
        clears the graph
        """
        self.node_cache.clear()
    
    def clear_touch_settings(self):
        """
        clear the touch settings for each of the nodes
        """

        for node in self.get_all_nodes():
            node.clear_touch()
    


    def add_edges(self, from_nodes_list, to_nodes_list):

        for from_node in from_nodes_list:
            for to_node in to_nodes_list:
                from_node.add_next_node(to_node)
                to_node.add_prev_node(from_node)

    def prune_edges(self, from_nodes_list, to_nodes_list):

        for from_node in from_nodes_list:
            for to_node in to_nodes_list:
                from_node.remove_next_node(to_node)
                to_node.remove_prev_node(from_node)
    

    def prune_node(self, node):
        self.prune_edges(node.get_prev_nodes(), list(node))
        self.prune_edges(list(node), node.get_next_nodes())
        self.node_cache.pop(node.get_loc_id())
    
        
    def get_gene_id(self):
        return self.gene_id


    def draw_graph(self, filename):
        
        ofh = open(filename, 'w')

        ofh.write("digraph G {\n")

        for node_id in self.node_cache:
            node = self.node_cache[node_id]
            node_seq = node.get_seq()
            gene_node_id = node.get_gene_node_id()
            next_nodes = node.get_next_nodes()

            ofh.write("{} [label=\"{}:Len{}:{}\"]\n".format(node.get_id(), gene_node_id, len(node_seq), node_seq))
            
            for next_node in next_nodes:
                ofh.write("{}->{}\n".format(node.get_id(), next_node.get_id()))

        ofh.write("}\n") # close it

        ofh.close()
    
