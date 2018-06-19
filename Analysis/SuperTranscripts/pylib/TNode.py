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

import TGraph

logger = logging.getLogger(__name__)


class TNode:
    """
    generic Trinity graph node object representing a node in the Trinity isoform reconstruction graph

    Node's are objects within a gene and can be shared among transcript isoforms.

    instance members include:

        tgraph : (TGraph obj) graph for the Trinity gene, which will hold the nodes.

        transcripts: list(str) names of the isoforms that contains this node.

        loc_node_id : (int) identifier of the node

        seq : (str)  nucleotide sequence for this node in the transcript

        len : (int)  length of the node sequence

        prev : (set) node objects connected as parental nodes in the graph

        next : (set) node objects connected as descendant nodes in the graph

    class members include:

        merged_nodeset_counter : (int) tracking nodes that get merged under squeeze operations.
        
    """
    
    node_cache = dict()

    merged_nodeset_counter = 0

    all_nodes_counter = 0


    def __init__(self, tgraph, transcript_id, loc_node_id, node_seq):
        """
        constructor, but don't use directly.... instead, use TGraph.get_node() factory function
        """
        
        if len(node_seq) == 0:
            raise RuntimeError("Error, TNode instantiation requires node sequence of length > 0")

        self.tgraph = tgraph
        self.transcripts = set()
        self.add_transcripts(transcript_id)
        self.loc_node_id = loc_node_id
        self.seq = node_seq
        self.len = len(node_seq)

        TNode.all_nodes_counter += 1
        self._id = TNode.all_nodes_counter
        
        #logger.info("{}\t{}".format(loc_node_id, node_seq))

        self.prev = set()
        self.next = set()
        self.stashed_prev = set() # for manipulation during topological sorting
        self.stashed_next = set() 

        self.touched = 0  
        self.dead = False
        self.topological_order = -1  # updated on topological sorting

    #########################
    ## various Node ID values
    #########################


    def get_id(self):
        # a private unique identifier for all nodes
        return self._id
        
    def get_loc_id(self):
        return self.loc_node_id

    def set_loc_id(self, loc_node_id):
        self.loc_node_id = loc_node_id
        
    
    def get_gene_id(self):
        return self.tgraph.get_gene_id()

    def get_gene_node_id(self):

        gene_id = self.get_gene_id()
        loc_id = self.get_loc_id()

        node_id = gene_id + "::" + loc_id
                
        return node_id

    def get_touched_val(self):
        return self.touched

    def is_dead(self):
        return(self.dead)

    def is_ancestral(self, node, visited=None):

        if visited is None:
            visited = set() #init round
        
        #logger.debug("is_ancestral search from {} of node {}".format(self, node))
        if node == self:
            #logger.debug("node is self")
            return True
        
        if node in self.prev:
            #logger.debug("node in self.prev")
            return True
        else:
            #logger.debug("continuing search")
            visited.add(self)
            #logger.debug("visited: {}".format(visited))
            for prev_node in self.prev:
                #logger.debug("cascading towards prev_node: {}".format(prev_node))
                if prev_node in visited:
                    #logger.debug("prev_node in visited")
                    pass
                else:
                    #logger.debug("prev_node not in visited")
                    found = prev_node.is_ancestral(node, visited)
                    if found:
                        return True
        return False


    def is_descendant(self, node, visited=None):

        if visited == None:
            visited = set()  # init round
        
        if node == self:
            return True

        if node in self.next:
            return True
        else:
            visited.add(self)
            for next_node in self.next:
                if next_node not in visited:
                    found = next_node.is_descendant(node, visited)
                    if found:
                        return True
        return False
    
    ## Other accessors

    def get_graph(self):
        return self.tgraph

    def get_seq(self):
        return self.seq


    def set_seq(self, seq):
        self.seq = seq


    def get_topological_order(self):
        return self.topological_order

    def set_topological_order(self, topo_order):
        self.topological_order = topo_order
    
    def get_transcripts(self):
        return self.transcripts

    def add_transcripts(self, transcript_name_or_set):
        if type(transcript_name_or_set) is set:
            self.transcripts.update(transcript_name_or_set)
        elif type(transcript_name_or_set) is str:
            self.transcripts.add(transcript_name_or_set)
        else:
            raise RuntimeError("Error, parameter must be a string or a set ")
        
    def get_prev_nodes(self):
        return set(self.prev)

    def get_next_nodes(self):
        return set(self.next)
    
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


    ## Touching nodes 
    def touch(self):
        self.touched += 1

    def untouch(self):
        self.touched -= 1

    def clear_touch(self):
        self.touched = 0


    def toString(self):
        txt = str("prev: " + str(self.get_prev_node_loc_ids()) +
                  ", me: " + str(self.get_loc_id()) +
                  ", next: " + str(self.get_next_node_loc_ids()) +
                  ", transcripts: " + str(self.transcripts) +
                  ", " + self.get_seq())

        if self.topological_order >= 0:
            txt += ", topo_order={}".format(self.topological_order)
        
        if self.dead:
            txt += " ** dead ** "
        
        return txt
    

    @classmethod
    def merge_nodes(cls, node_list):
        """
        Merges linear stretches of nodes into a single new node that has
        concatenated sequences of the input nodes
        """

        logger.debug("Merging nodes: {}".format(node_list))
        
        merged_node_seq = ""
        TNode.merged_nodeset_counter += 1
        merged_loc_node_id = "M{}".format(TNode.merged_nodeset_counter)

        # transcript list should be the intersection from nodes being merged (not the union)
        # because repeat nodes could be part of the merge.
        transcripts = node_list[0].get_transcripts()


        for node_obj in node_list:
            logger.debug("node being merge: {}".format(node_obj.toString()))
            seq = node_obj.get_seq()
            merged_node_seq += seq
            transcripts = transcripts.intersection(node_obj.get_transcripts())

        tgraph = node_list[0].get_graph()

        merged_node = TNode(tgraph, transcripts, merged_loc_node_id, merged_node_seq)
        
        return merged_node


    def is_burr(self):
        """
        returns true if node (x) is in this graphical context:

          X                               X
            \           or               /
         C-- A--?                   ?-- A--B

         where X dangles.


         So, X has only one parent or child and not otherwise connected in the graph.

         """

        if self.get_prev_nodes() and self.get_next_nodes():
            return False

        if len(self.get_prev_nodes()) > 1 or len(self.get_next_nodes()) > 1:
            return False

        # illustration above on left side
        if (len(self.get_next_nodes()) == 1
            and
            len(self.get_prev_nodes()) == 0
            and
            len(self.get_next_nodes().pop().get_prev_nodes()) > 1):

            return True

        # illustration above on right side
        if (len(self.get_next_nodes()) == 0
            and
            len(self.get_prev_nodes()) == 1
            and
            len(self.get_prev_nodes().pop().get_next_nodes()) > 1):

            return True

        # more complex structure
        return False
            
    
