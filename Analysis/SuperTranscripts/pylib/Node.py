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

    node_cache = dict()

    merged_nodeset_counter = 0

    def __init__(self, transcript_name, loc_node_id, node_seq):

        if len(node_seq) == 0:
            raise RuntimeError("Error, Node instantiation requires node sequence of length > 0")

        self.transcript_name = transcript_name
        self.loc_node_id = loc_node_id
        self.seq = node_seq
        self.len = len(node_seq)

        #logger.info("{}\t{}".format(loc_node_id, node_seq))

        self.prev = set()
        self.next = set()
        self.stashed_prev = set() # for manipulation during topological sorting
        self.stashed_next = set() 

    @classmethod
    def get_node(self, transcript_name, loc_node_id, node_seq):
        
        logger.debug("{}\t{}".format(loc_node_id, node_seq))
        
        if len(node_seq) == 0:
            raise RuntimeError("Error, non-zero length node_seq required for parameter")
        
        gene_name = Node.get_gene_name(transcript_name)
        node_id = self.get_node_id(gene_name, loc_node_id)
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
            node_obj = Node(transcript_name, loc_node_id, node_seq)
            Node.node_cache[ node_id ] = node_obj
            return node_obj


    @staticmethod
    def clear_node_cache():
        Node.node_cache.clear()
    
    @staticmethod
    def get_node_id(gene_name, loc_node_id):
        
        node_id = "::".join([gene_name, loc_node_id])
        return node_id

    
    def get_loc_id(self):
        return self.loc_node_id
    
    def get_seq(self):
        return self.seq

    def get_transcript_name(self):
        return self.transcript_name

    @staticmethod
    def get_gene_name(transcript_name):

        if re.match("^\^\^TRIN", transcript_name):
            # using internally specified topologically sorted graph
            return transcript_name

        (gene_name, count) = re.subn("_i\d+$", "", transcript_name)
        if count != 1:
            errmsg = "Error, couldn't extract gene_id from transcript_id: {}".format(transcript_name)
            logger.critical(errmsg)
            raise RuntimeError(errmsg)
        return gene_name
    
        
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
        merged_node_seq = ""
        Node.merged_nodeset_counter += 1
        merged_loc_node_id = "M{}".format(Node.merged_nodeset_counter)

        for node_obj in node_list:
            seq = node_obj.get_seq()
            merged_node_seq += seq

        merged_node = Node("trinity", merged_loc_node_id, merged_node_seq)

        return merged_node


