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
import GraphCycleException

logger = logging.getLogger(__name__)


class Topological_sort:

    """
    Implementation of Topological Sorting - as per described on wikipedia

    # https://en.wikipedia.org/wiki/Topological_sorting

    L ← Empty list that will contain the sorted elements
    S ← Set of all nodes with no incoming edge
    while S is non-empty do
        remove a node n from S
        add n to tail of L
        for each node m with an edge e from n to m do
            remove edge e from the graph
            if m has no other incoming edges then
                insert m into S
    if graph has edges then
        return error (graph has at least one cycle)
    else
        return L (a topologically sorted order)
    """

    @staticmethod
    def topologically_sort(nodes_list):

        logger.debug("Nodes to topo sort: " + str(nodes_list))
        
        L = list()  # empty list to contain the sorted elements
        S = list()  # set of all nodes with no incoming edge

        # populate S with nodes lacking incoming edges
        for node in nodes_list:
            if len(node.get_prev_nodes()) == 0:
                S.append(node)

        while len(S) > 0:
            node_n = S.pop(0)
            L.append(node_n)
            n_children = list(node_n.get_next_nodes()) # make a copy instead of using the reference directly.
            for n_child in n_children:
                n_child.stash_prev_node(node_n)
                node_n.stash_next_node(n_child)

                if len(n_child.get_prev_nodes()) == 0:
                    S.append(n_child)

        # see if there are cycles
        for node in nodes_list:
            if len(node.get_prev_nodes()) > 0 or len(node.get_next_nodes()) > 0:
                raise GraphCycleException("Graph has cycles!  offending node: " + node.toString())

        for node in L:
            node.restore_stashed_nodes()

        return L
    
