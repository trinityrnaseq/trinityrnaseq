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

class Node_alignment:

    """
    Object has two members:

        transcript_names = [ transA,
                             transB,
                             transC,
                             ...
                             ]

        aligned_nodes = [ [transA_node_1, transA_node_2, ... ],
                          [transB_node_1, transB_node_2, ... ],
                          [ None,         transC_node_1, ... ],  
                        ]

    Note, can have None at node positions to include gaps.

    """

    GAP = None

    def __init__(self, transcript_name_list, node_obj_matrix):
        self.transcript_names = transcript_name_list
        self.aligned_nodes = node_obj_matrix

    def get_transcript_names(self):
        # accessor
        return self.transcript_names

    def get_aligned_nodes(self):
        # accessor
        return self.aligned_nodes

    @staticmethod
    def get_single_seq_node_alignment(path_obj):
        """
        Factory method:
           constructs a Node_alignment object from a Node_path object

           mostly just reshaping the info for use with the multiple alignment methods.
        
        """
        
        node_list = list()
        for node_obj in path_obj.get_path():
            node_list.append(node_obj)

        self = Node_alignment([path_obj.get_transcript_name()], [node_list])

        return self


    @staticmethod
    def compute_number_common_nodes(align_A, align_B):
        """
        given to Node_alignment objects, counts the number of shared nodes
        """
        
        node_set_a = Node_alignment.get_node_set(align_A)
        node_set_b = Node_alignment.get_node_set(align_B)

        node_set_a = Node_alignment.get_node_loc_ids(node_set_a)
        node_set_b = Node_alignment.get_node_loc_ids(node_set_b)
                
        common_nodes = set.intersection(node_set_a, node_set_b)

        return common_nodes


    @staticmethod
    def get_node_loc_ids(node_set):
        """
        private static method
        gets the list of loc_id among all nodes in the set
        """
        
        loc_ids_set = set()
        for node in node_set:
            loc_id = node.get_loc_id()
            loc_ids_set.add(loc_id)

        return loc_ids_set
    

    @staticmethod
    def get_node_set(align_obj):
        """
        extracts a list of unique Node objects from the Node_alignment object
        """
        
        num_trans = len(align_obj)
        alignment_width = align_obj.width()

        node_set = set()
        
        for align_num in range(0,num_trans):
            for align_pos in range(0,alignment_width):
                node_obj = align_obj.aligned_nodes[ align_num ][ align_pos ]
                if node_obj is not None:
                    node_set.add(node_obj)

        return node_set


    def get_node_set_at_column_pos(self, col_pos):
        """
        At a given column of the Node_alignment, extracts the list of unique nodes
        """

        # FIXME: since we're not dealing with mismatched nodes, there really should only be one node here
        # that's shared among the different alignments
        # Need refactoring across the next few methods as well for the same reason.
        
        node_objs = set()
        for i in range(0, len(self)):
            node_obj = self.aligned_nodes[ i ][ col_pos ]
            if node_obj is not None:
                node_objs.add(node_obj)
        
        return node_objs

    def get_representative_column_node(self, col_pos):

        node_list = list(self.get_node_set_at_column_pos(col_pos))

        return node_list[ 0 ]
    

    def get_node_LIST_at_column_pos(self, col_pos):

        node_objs = list()
        for i in range(0, len(self)):
            node_obj = self.aligned_nodes[ i ][ col_pos ]
            node_objs.append(node_obj)
        
        return node_objs

    def get_node_occupancy_at_column_pos(self, col_pos):

        node_list = self.get_node_LIST_at_column_pos(col_pos)

        occupancy_list = list()
        for node in node_list:
            if node is None:
                occupancy_list.append(False)
            else:
                occupancy_list.append(True)

        return occupancy_list


    def append_node_to_each_entry(self, node_obj):

        for i in range(0, len(self)):
            self.aligned_nodes[ i ].append(node_obj)

    def append_node_according_to_occupancy_pattern(self, node_obj, occupancy_pattern):

        for i in range(0, len(self)):
            if occupancy_pattern[i] is True:
                self.aligned_nodes[ i ].append(node_obj)
            else:
                self.aligned_nodes[ i ].append(None)

        


    def add_column(self, column_node_list):
        num_alignments = len(self)
        if len(column_node_list) != num_alignments:
            raise RuntimeError("Error, column size differs from num_alignments")

        for i in range(0,num_alignments):
            self.aligned_nodes[ i ].append(column_node_list[ i ])
        
                    
    def __len__(self):
        """
        number of transcripts represented in the alignment
        """
        
        return(len(self.transcript_names))

    def width (self):
        """
        width of the alignment (number of columns)
        """
        return(len(self.aligned_nodes[0])) 

    
    def __repr__(self):

        num_transcripts = len(self.transcript_names)
        ret_text = "\n# Alignment obj contains: {} transcripts: {}\n\n".format(num_transcripts, ",".join(self.transcript_names))

        align_width = self.width()

        NODES_PER_LINE = 10

        # each alignment block
        for i in range(0, align_width, NODES_PER_LINE):

            # report alignment for each entry
            for j in range(0,num_transcripts):
                transcript_name = self.transcript_names[ j ]
                aligned_nodes_entry = self.aligned_nodes[ j ]

                ret_text += "{}".format(transcript_name)
                for x in range(i, i+NODES_PER_LINE):
                    if x >= align_width:
                        break
                    
                    ret_text += "\t{}".format(aligned_nodes_entry[ x ])

                ret_text += "\n" # end of current line

            ret_text += "\n" # spacer between alignment blocks
                        
            #ret_text += "Align [{}] trans {} : path {}".format(i, transcript_name, str(aligned_nodes_entry)) + "\n"

        return ret_text
    

    def squeeze(self):
        """
        merge unbranched nodes into single nodes
        """
        
        num_transcripts = len(self)
        width = self.width()

        node_obj_matrix = list()
        for i in range(0,num_transcripts):
            node_obj_matrix.append([])

        squeezed_alignment = Node_alignment(self.get_transcript_names(), node_obj_matrix)

        # walk the node list and merge unbranched stretches into single nodes
        block_breakpoints = []
        prev_col_node_set = self.get_node_occupancy_at_column_pos(0)
        for i in range(1,width):
            node_column_set = self.get_node_occupancy_at_column_pos(i)

            #print("Comparing {} to {} == {}".format(prev_col_node_set, node_column_set, prev_col_node_set == node_column_set))

            if node_column_set != prev_col_node_set:
                block_breakpoints.append(i)
            prev_col_node_set = node_column_set

        block_breakpoints.append(width)

        logger.debug("Block_breakpoints: {}".format(block_breakpoints))

        blocked_nodes = list()
        for i in range(0, width+1):
            if i in block_breakpoints:
                # found block terminator
                node_to_add = None
                if len(blocked_nodes) > 1:
                    node_to_add = TNode.TNode.merge_nodes(blocked_nodes)
                else:
                    node_to_add = blocked_nodes[ 0 ]

                blocked_node_occupancy = self.get_node_occupancy_at_column_pos(i-1)
                squeezed_alignment.append_node_according_to_occupancy_pattern(node_to_add, blocked_node_occupancy)

                blocked_nodes = list() # reinit
            
            # add to running block
            if i < width:
                blocked_nodes.append(self.get_representative_column_node(i))
        
        return squeezed_alignment


    def to_gene_fasta_and_gtf(self, gene_name):

        transcript_names = self.get_transcript_names()
        
        gene_seq = ""

        # init transcript gtf records
        transcript_to_gtf_lines = dict()

        transcript_to_malign = dict()

        for transcript_name in transcript_names:
            transcript_to_gtf_lines[ transcript_name ] = ""
            transcript_to_malign[ transcript_name ] = ""

        for i in range(0,self.width()):
            node_obj = self.get_representative_column_node(i)

            node_seq = node_obj.get_seq()
            if len(node_seq) == 0:
                raise RuntimeError("Error, node seq of length zero: node=" + str(node_obj))

            node_occupancy = self.get_node_occupancy_at_column_pos(i)

            pos_start = len(gene_seq) + 1
            gene_seq += node_obj.get_seq()
            pos_end = len(gene_seq)

            # include gtf record for transcripts
            for j in range(0,len(transcript_names)):
                transcript_name = transcript_names[ j ]
                if node_occupancy[ j ] is True:
                    # make gtf record
                    transcript_to_gtf_lines[ transcript_name ] += "\t".join([gene_name, "Trinity_gene", "exon",
                                                                            str(pos_start), str(pos_end), '.', '+', '.',
                                                                            "gene_id \"{}\"; transcript_id \"{}\"\n".format(
                                                                                gene_name, transcript_name) ] )
                    transcript_to_malign[ transcript_name ] += node_seq
                else:
                    for x in range(0,len(node_seq)):
                        transcript_to_malign[ transcript_name ] += '.'
                
        
        # build mini-gtf section
        gene_gtf = "\n".join(transcript_to_gtf_lines.values())

        return (gene_seq, gene_gtf, transcript_to_malign)
        
