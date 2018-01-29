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
import Trinity_util

logger = logging.getLogger(__name__)


class Node_path:
    """
    Object representation of the connected set of Node objects that represent the reconstructed isoforms graph traversal

    Instance members:

        transcript_name : (str)  name of the isoform

        node_obj_list : (list) of Node objects

    """

    
    def __init__(self, tgraph, transcript_name, path_string, sequence):
        """
        constructor, instantiates Node_path and builds vertices in the graph
        """
        
        self.transcript_name = transcript_name
        self.node_obj_list = list()

        node_descr_list = re.findall("\d+:\d+\-\d+", path_string)

        obj_node_list = list()
        for node_descr in node_descr_list:
            (loc_node_id, node_coord_range) = node_descr.split(":")
            (lend,rend) = node_coord_range.split("-")
            lend = int(lend)
            rend = int(rend)

            # use factory call to instantiate node objects:
            node_obj = tgraph.get_node(transcript_name,
                                     loc_node_id, sequence[lend:rend+1]) # coords in path were already zero-based
            
            self.node_obj_list.append(node_obj)

    
    def get_transcript_name(self):
        return self.transcript_name

    def get_path(self):
        return self.node_obj_list
    

    def __repr__(self):
        node_str_list = list()
        for node in self.node_obj_list:
            node_str_list.append(str(node))

        path_str = "--".join(node_str_list)

        return path_str
        
