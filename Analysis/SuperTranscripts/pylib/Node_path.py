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

        # 1st node is special (full kmer prefix included)
        first_kmer_flag = False
        obj_node_list = list()
        for node_descr in node_descr_list:
            (loc_node_id, node_coord_range) = node_descr.split(":")
            (lend,rend) = node_coord_range.split("-")
            lend = int(lend)
            rend = int(rend)

            if not first_kmer_flag:
                first_kmer_flag = True
                loc_node_id += 'fst'
            
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
        

    @staticmethod
    def adjust_for_fst_nodes(tgraph, node_path_list):
        """
        fst nodes will have an extra 5' sequence as compared to the corresponding non-fst nodes.

        If both the fst and non-fst version of the node exist, must modify the fst nodes so that
        they are separated from their 5' extension, and the core of the node (suffix) is shared.

        input: TGraph obj, list of node_path objects.

        The node_path objects are modified in-place as needed.
        A fst-node will be truncated to the unique prefix and the non-fst node will be integrated into the path.

        returns the node_path_list with any required adjustments

        """

        # get list of fst nodes requiring adjustment

        fst_nodes_require_adj = list()

        nodes = tgraph.get_all_nodes()
        for node in nodes:
            node_id = node.get_loc_id()
            if re.search("fst", node_id):
                core_node_id = re.sub("fst", "", node_id)
                core_node = tgraph.retrieve_node(core_node_id)
                if core_node is not None:
                    fst_nodes_require_adj.append( (node, core_node) )

        if not fst_nodes_require_adj:
            # nothing to do
            logger.debug("no FST nodes to adjust")
            return node_path_list


        logger.debug("Adjusting FST nodes: {}".format(fst_nodes_require_adj))

        old_fst_node_to_new_fst_nodes = dict()
        fst_nodes_to_delete = list()

        # perform node modifications:
        for (fst_node, core_node) in fst_nodes_require_adj:

            fst_node_seq = fst_node.get_seq()
            core_node_seq = core_node.get_seq()

            prefix_endpt = fst_node_seq.index(core_node_seq)



            if prefix_endpt == 0:
                assert fst_node_seq == core_node_seq, "Error, prefix starts at first position but sequences are not equivalent"
                core_node.add_transcripts(fst_node.get_transcripts())
                old_fst_node_to_new_fst_nodes[fst_node] = [core_node]
                fst_nodes_to_delete.append(fst_node)

            else:
                prefix_string = fst_node_seq[0:prefix_endpt]
                logger.debug("FST-SEQ-EXTRACTION\n\nFSTseq:\n{}\n\nCOREseq:\n{}\n\nPREFIXseq:\n{}\n\n".format(fst_node_seq, core_node_seq, prefix_string))
                core_node.add_transcripts(fst_node.get_transcripts())
                old_fst_node_to_new_fst_nodes[fst_node] = [fst_node, core_node]
                fst_node.set_seq(prefix_string)

        # now perform node path updates
        for node_path in node_path_list:
            nodes = node_path.get_path()

            first_node = nodes[0]
            if first_node in old_fst_node_to_new_fst_nodes:
                # must replace
                replacement_node_list = old_fst_node_to_new_fst_nodes[first_node]
                if len(replacement_node_list) == 1:
                    # swap it out
                    nodes[0] = replacement_node_list[0]
                elif len(replacement_node_list) == 2:
                    nodes[0] = replacement_node_list[1]
                    nodes.insert(0, replacement_node_list[0])
                else:
                    raise RuntimeError("shouldn't get here")
                
        # purge the nodes targeted for deletion
        for node in fst_nodes_to_delete:
            tgraph.prune_node(node)

        return node_path_list

        
