#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

import TGraph
import TNode
import Node_alignment
from Topological_sort import Topological_sort
import TGLOBALS

logger = logging.getLogger(__name__)



MAX_MM_RATE = 0.05 


class Compact_graph_whole:

    def __init__(self):

        pass

    
    def compact_unbranched(self, tgraph):

        all_nodes = Topological_sort.topologically_sort(tgraph.get_all_nodes())
        
        for node in all_nodes:

            if node.is_dead():
                continue

            logger.debug("compact_unbranched( {} )".format(node.toString()))
            prev_nodes = sorted(node.get_prev_nodes())
            if len(prev_nodes) == 1:
                # merge them
                prev_node = prev_nodes[0]
                if len(prev_node.get_next_nodes()) != 1:
                    # prev node is downward branched... no dice
                    logger.debug("prev node: {} id downward branched, so can't compact. {}".format(prev_node, prev_node.toString()))
                    continue

                logger.debug("linear compaction of nodes {} and {}".format(prev_node, node))
                node.add_transcripts(prev_node.get_transcripts())
                node.set_seq(prev_node.get_seq() + node.get_seq())
                tgraph.add_edges(prev_node.get_prev_nodes(), [node])
                tgraph.prune_node(prev_node)


    def compact_graph(self, tgraph, num_allowed_variants):

        compacted_flag = True

        round = 0
        while compacted_flag:
            
            round += 1
            logger.debug("\n\t### Num allowed variants: {}, Round: {}".format(num_allowed_variants, round))

            tgraph.clear_touch_settings() # start fresh
            compacted_flag = False # reset for this round

            sorted_nodes = Topological_sort.topologically_sort(tgraph.get_all_nodes())

            if TGLOBALS.DEBUG:
                dot_filename = "ladeda.compacted_whole.BEGIN.num_allowed_{}.round_{}.dot".format(num_allowed_variants, round)
                tgraph.draw_graph(dot_filename)
                logger.debug("-wrote dot: {}".format(dot_filename))
            
            for node in sorted_nodes:

                if node.get_touched_val() > 0:
                    continue

                logger.debug("\n\n///////// R[{}]  Exploring compaction at node: {}\n\n".format(round, node.toString()))
                
                # try compact upward
                prev_nodes = node.get_prev_nodes()
                if len(prev_nodes) > 1 and self.untouched(prev_nodes) and self.all_have_lower_topological_orderings(node, prev_nodes):
                    if self.compact_upward(prev_nodes, num_allowed_variants):
                        compacted_flag = True

                next_nodes = node.get_next_nodes()
                if len(next_nodes) > 1 and self.untouched(next_nodes) and self.all_have_higher_topological_orderings(node, next_nodes):
                    if self.compact_downward(next_nodes, num_allowed_variants):
                        compacted_flag = True

            if compacted_flag:
                if TGLOBALS.DEBUG:
                    dot_filename = "ladeda.compacted_whole.END.num_allowed_{}.round_{}.dot".format(num_allowed_variants, round)
                    tgraph.draw_graph(dot_filename)
                    logger.debug("-wrote dot: {}".format(dot_filename))

                self.compact_unbranched(tgraph)


    def untouched(self, node_list):
        for node in node_list:
            if node.get_touched_val() > 1:
                return False

        return True


    def all_have_lower_topological_orderings(self, node, other_nodes_list):
        topo_order = node.get_topological_order()
        if topo_order < 0:
            raise RuntimeError("Error, must run topological ordering first")

        for other_node in other_nodes_list:
            if other_node.get_topological_order() > topo_order:
                return False
        return True


    def all_have_higher_topological_orderings(self, node, other_nodes_list):
        topo_order = node.get_topological_order()
        if topo_order < 0:
            raise RuntimeError("Error, must run topological ordering first")

        for other_node in other_nodes_list:
            if other_node.get_topological_order() < topo_order:
                return False
        return True
    

    ####################
    ## Upward Compaction
    
    def compact_upward(self, prev_nodes, num_allowed_variants):

        prev_nodes = sorted(prev_nodes)

        logger.debug("compact_upward: {} max_allowed_variants: {}".format(prev_nodes, num_allowed_variants))

        for i in range(0, len(prev_nodes)-1):
            for j in range(i+1, len(prev_nodes)):
                if self.compact_upward_node_pair(prev_nodes[i], prev_nodes[j], num_allowed_variants):
                    return True

        return False

    def compact_upward_node_pair(self, node_A, node_B, num_allowed_variants):
        logger.debug("compact_upward_node_pair({}, {}, num_allowed_variants: {}".format(node_A, node_B, num_allowed_variants))
        
        # ensure they're on parallel paths in the graph.
        if (node_A.is_ancestral(node_B) or node_B.is_ancestral(node_A)):
            logger.debug("-not parallel paths. excluding compaction of {} and {}".format(node_A, node_B))
            return False
        logger.debug("-not ancestral: {}, {}".format(node_A, node_B))
        
        # switch A,B so A is the shorter sequence node

        if len(node_B.get_seq()) < len(node_A.get_seq()):
            (node_A, node_B) = (node_B, node_A)

        seqA = node_A.get_seq()
        seqB = node_B.get_seq()

        shorter_len = len(seqA)

        # reverse the strings and check number of mismatches
        seqA_rev = seqA[::-1]
        seqB_rev = seqB[::-1]

        num_mm = 0
        for i in range(0, shorter_len):
            if seqA_rev[i] != seqB_rev[i]:
                num_mm += 1

        if num_mm > num_allowed_variants:
            logger.debug("num_mm: {} exceeds {} for node pair: {} and {}".format(num_mm, num_allowed_variants, node_A, node_B))
            return False

        if num_mm / shorter_len > MAX_MM_RATE:
            logger.debug("num_mm: {} / len: {} exceeds max mm rate: {} for node pair: {} and {}".format(num_mm, shorter_len,
                                                                                                        MAX_MM_RATE, node_A, node_B))
            return False


        # go ahead and merge the two in their region of common overlap
        tgraph = node_A.get_graph()

        ##### operations needed:
        #
        #         A
        #           \
        #            -- X
        #           / 
        #         B
        #
        #   -add B transcripts to A
        #   -adjust all B-downward edges to A
        #   -if B is longer than A, set upward to A and trim B part to unique region
        #      otherwise, remove B altogether.
        #

        logger.debug("node_A before mods: {}".format(node_A.toString()))
        logger.debug("node_B before mods: {}".format(node_B.toString()))

        node_A.add_transcripts(node_B.get_transcripts())
        tgraph.add_edges( [node_A], node_B.get_next_nodes() )
        tgraph.prune_edges( [node_B], node_B.get_next_nodes() )

        node_A.touch()
        node_B.touch()


        if len(seqB) > shorter_len:
            delta = len(seqB) - shorter_len
            seqB = seqB[0:delta]
            node_B.set_seq(seqB)
            tgraph.add_edges([node_B], [node_A])
        else:
            # remove node_B prev edges and attach them to A
            tgraph.add_edges(node_B.get_prev_nodes(), [node_A])
            # remove B altogether
            tgraph.prune_node(node_B)

        logger.debug("node_A after mods: {}".format(node_A.toString()))
        logger.debug("node_B after mods: {}".format(node_B.toString()))

        logger.debug("\n\n\t*** compacted nodes: {} and {} ***".format(node_A, node_B))

        return True

    ######################
    ## Downward compaction

    def compact_downward(self, next_nodes, num_allowed_variants):

        next_nodes = sorted(next_nodes)
        logger.debug("compact_downward: {} max_allowed_variants: {}".format(next_nodes, num_allowed_variants))

        for i in range(0, len(next_nodes)-1):
            for j in range(i+1, len(next_nodes)):
                if self.compact_downward_node_pair(next_nodes[i], next_nodes[j], num_allowed_variants):
                    return True

        return False

    def compact_downward_node_pair(self, node_A, node_B, num_allowed_variants):

        # ensure they're on parallel paths in the graph.
        if (node_A.is_descendant(node_B) or node_B.is_descendant(node_A)):
            logger.debug("-not parallel paths. excluding compaction of {} and {}".format(node_A, node_B))
            return False

        # switch A,B so A is the shorter sequence node

        if len(node_B.get_seq()) < len(node_A.get_seq()):
            (node_A, node_B) = (node_B, node_A)

        seqA = node_A.get_seq()
        seqB = node_B.get_seq()

        shorter_len = len(seqA)


        num_mm = 0
        for i in range(0, shorter_len):
            if seqA[i] != seqB[i]:
                num_mm += 1

        if num_mm > num_allowed_variants:
            logger.debug("num_mm: {} exceeds {} for node pair: {} and {}".format(num_mm, num_allowed_variants, node_A, node_B))
            return False

        if num_mm / shorter_len > MAX_MM_RATE:
            logger.debug("num_mm: {} / len: {} exceeds max mm rate: {} for node pair: {} and {}".format(num_mm, shorter_len,
                                                                                                        MAX_MM_RATE, node_A, node_B))
            return False


        # go ahead and merge the two in their region of common overlap
        tgraph = node_A.get_graph()

        ##### operations needed:
        #
        #                   A
        #                 /
        #            -- X
        #                 \ 
        #                   B
        #
        #   -add B transcripts to A
        #   -adjust all B-downward edges to A
        #   -if B is longer than A, set downward to A and trim B part to unique region
        #      otherwise, remove B altogether.
        #

        logger.debug("node_A before mods: {}".format(node_A.toString()))
        logger.debug("node_B before mods: {}".format(node_B.toString()))

        node_A.add_transcripts(node_B.get_transcripts())
        tgraph.add_edges(node_B.get_prev_nodes(), [node_A])
        tgraph.prune_edges(node_B.get_prev_nodes(), [node_B])

        node_A.touch()
        node_B.touch()

        if len(seqB) > shorter_len:
            seqB = seqB[shorter_len:]
            node_B.set_seq(seqB)
            tgraph.add_edges([node_A], [node_B])
        else:
            # remove node_B prev edges and attach them to A
            tgraph.add_edges([node_A], node_B.get_next_nodes())
            # remove B altogether
            tgraph.prune_node(node_B)

        logger.debug("node_A after mods: {}".format(node_A.toString()))
        logger.debug("node_B after mods: {}".format(node_B.toString()))

        logger.debug("\n\n\t*** compacted nodes: {} and {} ***".format(node_A, node_B))

        return True


