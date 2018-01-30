#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, sys, re
import logging

import TGraph
import TNode
import Node_alignment
import Compact_graph_whole


logger = logging.getLogger(__name__)
#logger.addHandler(logging.NullHandler())


class Compact_graph_partial(Compact_graph_whole.Compact_graph_whole):

    def __init__(self):

        pass

    
    def compact_graph(self, tgraph):

        compacted_flag = True

        round = 0
        while compacted_flag:

            round += 1
            logger.debug("\n\t### Partial graph compaction, Round: {}".format(round))

            tgraph.clear_touch_settings() # start fresh
            compacted_flag = False # reset for this round

            for node in tgraph.get_all_nodes():

                if node.get_touched_val() > 0:
                    continue

                # try compact upward
                prev_nodes = node.get_prev_nodes()
                if len(prev_nodes) > 1 and self.untouched(prev_nodes):
                    if self.compact_upward(prev_nodes):
                        compacted_flag = True

                next_nodes = node.get_next_nodes()
                if len(next_nodes) > 1 and self.untouched(next_nodes):
                    if self.compact_downward(next_nodes):
                        compacted_flag = True

            if compacted_flag:
                dot_filename = "ladeda.compacted_partial.round_{}.dot".format(round)
                tgraph.draw_graph(dot_filename)
                logger.debug("-wrote dot: {}".format(dot_filename))

                self.compact_unbranched(tgraph)


    ####################
    ## Upward Compaction
    
    def compact_upward(self, prev_nodes):

        logger.debug("compact_upward: Partial {}".format(prev_nodes))

        prev_nodes = list(prev_nodes)

        for i in range(0, len(prev_nodes)-1):
            for j in range(i+1, len(prev_nodes)):
                if self.compact_upward_node_pair(prev_nodes[i], prev_nodes[j]):
                    return True

        return False

    def compact_upward_node_pair(self, node_A, node_B):

        logger.debug("compact_upward_node_pair({},{})".format(node_A, node_B))
        # ensure they're on parallel paths in the graph.
        if (node_A.is_ancestral(node_B) or node_B.is_ancestral(node_A)):
            logger.debug("-not parallel paths. excluding compaction of {} and {}".format(node_A, node_B))
            return False

        # switch A,B so A is the shorter sequence node
        
        if len(node_B.get_seq()) < len(node_A.get_seq()):
            (node_A, node_B) = (node_B, node_A)

        seqA = node_A.get_seq()
        seqB = node_B.get_seq()
        
        shorter_len = len(seqA)

        # reverse the strings and check number of mismatches
        seqA_rev = seqA[::-1]
        seqB_rev = seqB[::-1]

        num_matches = 0
        for i in range(0, shorter_len):
            if seqA_rev[i] == seqB_rev[i]:
                num_matches += 1
            else:
                break

        if num_matches == 0:
            logger.debug("no matching suffix between {} and {}".format(node_A, node_B))
            return False
        
        
        # go ahead and merge the two in their region of common overlap
        tgraph = node_A.get_graph()

        shared_seq = seqA[-1*num_matches:]
        uniqA_seq = seqA[0:(len(seqA)-num_matches)]
        uniqB_seq = seqB[0:(len(seqB)-num_matches)]

        logger.debug("Suffices:\nseqA:{}\nSeqB:{}\nShared:{}".format(seqA, seqB, shared_seq))
        

        assert len(seqA) == len(shared_seq) + len(uniqA_seq)
        assert len(seqB) == len(shared_seq) + len(uniqB_seq)
        
        ##### operations needed:
        #
        #         A                  A
        #           \                 \
        #             X       -->      C --X
        #           /                 /
        #         B                  B
        #
        #   - create C and set to partial shared sequence
        #   - add B and A transcripts to C
        #   - A and B link to C only
        #   - C recapitulates all outgoing edges from  A and B

        logger.debug("node_A before mods: {}".format(node_A.toString()))
        logger.debug("node_B before mods: {}".format(node_B.toString()))

        combined_transcripts = node_A.get_transcripts().union(node_B.get_transcripts())
        combined_loc_id = "Upart:" + node_A.get_loc_id() + "," + node_B.get_loc_id()
        node_C = tgraph.get_node(combined_transcripts, combined_loc_id, shared_seq)

        node_A.set_seq(uniqA_seq)
        node_B.set_seq(uniqB_seq)

        # reset next nodes from A,B to C
        all_next_nodes = node_A.get_next_nodes().union(node_B.get_next_nodes())
        tgraph.prune_edges([node_A], node_A.get_next_nodes())
        tgraph.prune_edges([node_B], node_B.get_next_nodes())
        tgraph.add_edges([node_C], all_next_nodes)
        node_C.touch()
        
        if len(uniqA_seq) == 0:
            # no longer need this now
            tgraph.add_edges(node_A.get_prev_nodes(), [node_C])
            tgraph.prune_node(node_A)
        else:
            tgraph.add_edges([node_A], [node_C])
            node_A.touch()

        if len(uniqB_seq) == 0:
            # no longer need it
            tgraph.add_edges(node_B.get_prev_nodes(), [node_C])
            tgraph.prune_node(node_B)
        else:
            tgraph.add_edges([node_B], [node_C])
            node_B.touch()
            
                
        logger.debug("node_A after mods: {}".format(node_A.toString()))
        logger.debug("node_B after mods: {}".format(node_B.toString()))
        logger.debug("node_C after mods: {}".format(node_C.toString()))
        
        logger.debug("\n\n\t*** partially UP-compacted nodes: {} and {} + {} ***".format(node_A, node_B, node_C))
        
        return True

    ######################
    ## Downward compaction

    def compact_downward(self, next_nodes):

        logger.debug("compact_downward: Partial {}".format(next_nodes))

        next_nodes = list(next_nodes)
        for i in range(0, len(next_nodes)-1):
            for j in range(i+1, len(next_nodes)):
                if self.compact_downward_node_pair(next_nodes[i], next_nodes[j]):
                    return True

        return False

    def compact_downward_node_pair(self, node_A, node_B):

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


        num_matches = 0
        for i in range(0, shorter_len):
            if seqA[i] == seqB[i]:
                num_matches += 1
            else:
                break


        if num_matches == 0:
            logger.debug("Error, no prefix match between {} and {}".format(node_A, node_B))
            return False

        
        # go ahead and merge the two in their region of common overlap
        tgraph = node_A.get_graph()

        shared_seq = seqA[0:num_matches]
        uniqA_seq = seqA[num_matches:]
        uniqB_seq = seqB[num_matches:]
        
        logger.debug("Prefixes:\nseqA:{}\nSeqB:{}\nShared:{}".format(seqA, seqB, shared_seq))

        ##### operations needed:
        #
        #                  A                          A
        #                 /                          /
        #            -- X           -->     -- X -- C
        #                 \                          \ 
        #                  B                          B
        #
        #  - find shared prefix of A,B
        #  - create C for shared prefix, shorten A,B
        #  - link A,B down from C
        #  - C takes on all prev nodes from A,B

        logger.debug("node_A before mods: {}".format(node_A.toString()))
        logger.debug("node_B before mods: {}".format(node_B.toString()))

        combined_transcripts = node_A.get_transcripts().union(node_B.get_transcripts())
        combined_loc_id = "Dpart:" + node_A.get_loc_id() + "," + node_B.get_loc_id()
        node_C = tgraph.get_node(combined_transcripts, combined_loc_id, shared_seq)

        node_A.set_seq(uniqA_seq)
        node_B.set_seq(uniqB_seq)

        # reset prev nodes from A,B to C
        all_prev_nodes = node_A.get_prev_nodes().union(node_B.get_prev_nodes())
        tgraph.prune_edges(node_A.get_prev_nodes(), [node_A])
        tgraph.prune_edges(node_B.get_prev_nodes(), [node_B])
        tgraph.add_edges(all_prev_nodes, [node_C])
        node_C.touch()

        
        if len(uniqA_seq) == 0:
            # no longer need this now
            tgraph.add_edges([node_C], node_A.get_next_nodes())
            tgraph.prune_node(node_A)
        else:
            tgraph.add_edges([node_C], [node_A])
            node_A.touch()

        if len(uniqB_seq) == 0:
            # no longer need it
            tgraph.add_edges([node_C], node_B.get_next_nodes())
            tgraph.prune_node(node_B)
        else:
            tgraph.add_edges([node_C], [node_B])
            node_B.touch()

        
        logger.debug("node_A after mods: {}".format(node_A.toString()))
        logger.debug("node_B after mods: {}".format(node_B.toString()))
        logger.debug("node_C after mods: {}".format(node_C.toString()))
        
        logger.debug("\n\n\t*** partially compacted nodes: {} and {} + {} ***".format(node_A, node_B, node_C))

        return True


