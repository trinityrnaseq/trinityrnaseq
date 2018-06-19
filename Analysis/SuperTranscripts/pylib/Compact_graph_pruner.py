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
    
    
    
