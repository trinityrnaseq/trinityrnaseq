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

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


class Node:
    """
    generic graph node object representing a node in the Trinity isoform reconstruction graph

    Node's are objects within a gene and can be shared among transcript isoforms.

    instance members include:

        gene_id : (str) name of the Trinity gene that the node corresponds to.

        loc_node_id : (int) identifier of the node

        seq : (str)  nucleotide sequence for this node in the transcript

        len : (int)  length of the node sequence

        prev : (set) node objects connected as parental nodes in the graph

        next : (set) node objects connected as descendant nodes in the graph

    class members include:

        node_cache : (dict) stores all nodes instantiated via the get_node() factory constructor.

        merged_nodeset_counter : (int) tracking nodes that get merged under squeeze operations.
        
    """
    
    node_cache = dict()

    merged_nodeset_counter = 0

    def __init__(self, gene_id, loc_node_id, node_seq):
        """
        constructor, but don't use directly.... instead, use get_node() factory function below
        """
        
        if len(node_seq) == 0:
            raise RuntimeError("Error, Node instantiation requires node sequence of length > 0")

        self.gene_id = gene_id
        self.loc_node_id = loc_node_id
        self.seq = node_seq
        self.len = len(node_seq)

        #logger.info("{}\t{}".format(loc_node_id, node_seq))

        self.prev = set()
        self.next = set()
        self.stashed_prev = set() # for manipulation during topological sorting
        self.stashed_next = set() 
        


    @classmethod
    def get_node(cls, gene_id, loc_node_id, node_seq):
        
        """
        Instantiates Node objects, and stores them in a graph.

        *** use this method for instantiating all Node objects ***

        use Node.clear_node_cache() to clear the graph
                
        """
                
        logger.debug("{}\t{}".format(loc_node_id, node_seq))
        
        if len(node_seq) == 0:
            raise RuntimeError("Error, non-zero length node_seq required for parameter")

        node_id = cls._construct_gene_node_id(gene_id, loc_node_id)
        if node_id in Node.node_cache:
            node_obj = Node.node_cache[ node_id ]
            if node_obj.seq != node_seq:
                ## FIXME: when internal node is found as a terminal node, we're using the k-1mer now for the terminal node too.
                logger.debug("WARNING: have conflicting node sequences for node_id: {}\n".format(node_id) +
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
            node_obj = Node(gene_id, loc_node_id, node_seq)
            Node.node_cache[ node_id ] = node_obj
            return node_obj


    @staticmethod
    def clear_node_cache():
        """
        clears the graph
        """
        Node.node_cache.clear()
    
    @staticmethod
    def _construct_gene_node_id(gene_id, loc_node_id):
        """
        builds a node identifier as a combination of the gene_name and loc_node_id
        """
        node_id = "::".join([gene_id, loc_node_id])
        return node_id

    
    def get_loc_id(self):
        return self.loc_node_id
    
    def get_seq(self):
        return self.seq

    def get_transcript_name(self):
        return self.transcript_name

    def get_prev_nodes(self):
        return self.prev

    def get_next_nodes(self):
        return self.next
    

    def add_next_node(self, next_node_obj):
        self.next.add(next_node_obj)

    def remove_next_node(self, remove_node_obj):
        self.next.remove(remove_node_obj)

    def stash_next_node(self, stash_node_obj):
        self.remove_next_node(stash_node_obj)
        self.stashed_next.add(stash_node_obj)

    def add_prev_node(self, prev_node_obj):
        self.prev.add(prev_node_obj)

    def remove_prev_node(self, remove_node_obj):
        self.prev.remove(remove_node_obj)

    def stash_prev_node(self, stash_node_obj):
        self.remove_prev_node(stash_node_obj)
        self.stashed_prev.add(stash_node_obj)


    def restore_stashed_nodes(self):
        self.prev.update(self.stashed_prev)
        self.stashed_prev = set()

        self.next.update(self.stashed_next)
        self.stashed_next = set()

    def get_prev_node_loc_ids(self):
        loc_ids = list()
        for node in self.get_prev_nodes():
            loc_ids.append(node.get_loc_id())
        return loc_ids

    def get_next_node_loc_ids(self):
        loc_ids = list()
        for node in self.get_next_nodes():
            loc_ids.append(node.get_loc_id())
        return loc_ids
    
    def __repr__(self):
        return(self.loc_node_id)


    def toString(self):
        txt = str("prev: " + str(self.get_prev_node_loc_ids()) +
                  ", me: " + str(self.get_loc_id()) +
                  ", next: " + str(self.get_next_node_loc_ids()))

        return txt
    
        
    @classmethod
    def merge_nodes(cls, node_list):
        """
        Merges linear stretches of nodes into a single new node that has
        concatenated sequences of the input nodes
        """
        
        merged_node_seq = ""
        Node.merged_nodeset_counter += 1
        merged_loc_node_id = "M{}".format(Node.merged_nodeset_counter)

        for node_obj in node_list:
            seq = node_obj.get_seq()
            merged_node_seq += seq

        merged_node = Node("trinity", merged_loc_node_id, merged_node_seq)

        return merged_node


