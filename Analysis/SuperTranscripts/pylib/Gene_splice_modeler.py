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
import TNode
import Node_path
import Node_alignment
from GraphCycleException import GraphCycleException 
import Topological_sort
import DP_matrix

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

class Gene_splice_modeler:

    """
    Builds supertranscipts.

    object instance members:

        gene_id : str

        alignments : list of Node_alignment objects

    """
    
    def __init__(self, gene_id, node_path_obj_list):

        """
        initialize alignments list with simple single 'alignment' objects with
        each path as an individual alignment with just its path nodes.

        params:

        gene_id : str

        node_path_obj_list : list of Node_path objects, each Node_path corresponding to an individual Trinity isoform

        """
        
        self.gene_id = gene_id
        self.alignments = list()

        logger.debug("Gene_splice_modeler inputs: {}".format(node_path_obj_list))
        
        for node_path_obj in node_path_obj_list:
            transcript_name = node_path_obj.get_transcript_name()
            alignment_obj = Node_alignment.Node_alignment.get_single_seq_node_alignment(node_path_obj)

            self.alignments.append(alignment_obj)

        

    def get_gene_id(self):
        return self.gene_id

    
    def build_splice_model(self):
        """
        method to construct the super transcript.

        Tries 2 approaches:
            a.  If there isn't an obvious repetitive node structure and so the graph formas a DAG,
                we build a splice graph and perform topological sorting of the nodes.
            b.  If there is some repetitive structure, we resort to performing a multiple alignment-based method to
                organize relationships among nodes in isoforms, and the multiple alignment produces the linear ordering
                for the supertranscript.

        """

        
        if not self.alignment_contains_repeat_node():
            # no obvious cycles
            try:
                return self.topological_order_splice_model()
            except GraphCycleException:
                # have a more complex cycle here...
                # try again w/ mult align approach.
                return self.multiple_alignment_splice_model()
            
        else:
            return self.multiple_alignment_splice_model()


    def alignment_contains_repeat_node(self):

        for alignment in self.alignments:
            loc_ids = set()
            for i in range(0, alignment.width()):
                node_obj = alignment.get_representative_column_node(i)
                loc_id = node_obj.get_loc_id()
                if loc_id in loc_ids:
                    return True
                loc_ids.add(loc_id)

        return False


    def topological_order_splice_model(self):
        """
        Build supertranscript using simpler topological sorting of the nodes.
        """
        
        logger.debug("\tusing topological sort method.\n");
        gene_id = self.get_gene_id()
        
        ## make a generic graph.
        graph = TGraph.TGraph(gene_id)
        for alignment in self.alignments:
            logger.debug("topological_order_splice_model, input alignment: " + str(alignment))
            node_list = alignment.get_aligned_nodes()[0] # should be unaligned here, so just ordered path list.
            transcript_name = alignment.get_transcript_names()[0]
            logger.debug("topological_order_splice_model, node list: " + str(node_list))
            for i in range(0, len(node_list)):
                node_obj = node_list[i]
                loc_id = node_obj.get_loc_id()
                generic_node = graph.get_node(transcript_name, loc_id, node_obj.get_seq()) # rely on Node class caching system
                logger.debug("generic node: " + str(generic_node))
                
                if i > 0:
                    # set prev node info
                    prev_node_obj = node_list[i-1]
                    prev_generic_node = graph.get_node(transcript_name, prev_node_obj.get_loc_id(), prev_node_obj.get_seq())
                    generic_node.add_prev_node(prev_generic_node)

                if i < len(node_list) - 1:
                    next_node_obj = node_list[i+1]
                    next_generic_node = graph.get_node(transcript_name, next_node_obj.get_loc_id(), next_node_obj.get_seq())
                    generic_node.add_next_node(next_generic_node)

        logger.debug("Before sorting nodes: " + str(graph))

        topologically_sorted_nodes = Topological_sort.Topological_sort.topologically_sort(graph.get_all_nodes())

        logger.debug("Topologically sorted nodes: " + str(topologically_sorted_nodes))
        
        # index loc node ids
        aligned_loc_id_pos = dict()
        for i in range(0, len(topologically_sorted_nodes)):
            loc_id = topologically_sorted_nodes[i].get_loc_id()
            aligned_loc_id_pos[loc_id] = i


        new_alignments = list()
        transcript_ids = list()
        for alignment in self.alignments:
            transcript_ids.append(alignment.get_transcript_names()[0]) # really should only be one here.
            new_alignment = [None for i in topologically_sorted_nodes]
            for node in alignment.get_aligned_nodes()[0]:
                loc_id = node.get_loc_id()
                new_idx = aligned_loc_id_pos[loc_id]
                new_alignment[new_idx] = node
            new_alignments.append(new_alignment)

        splice_graph_model = Node_alignment.Node_alignment(gene_id, transcript_ids, new_alignments)

        logger.debug("Splice graph model: " + str(splice_graph_model))

        return splice_graph_model
    
    def multiple_alignment_splice_model(self):
        """
        Multiple alignment algorithm for dealing with repeat nodes:
        For each best matching pair of transcripts (or aligned transcripts),
        perform alignment, and replace aligned pair with a single alignment object.
        """
        
        logger.debug("\tusing mult alignment method.\n");
                    
        alignments = self.alignments

        if len(alignments) == 1:
            # no alignment is necessary.
            return alignments[0]
        
        # determine initial path similarity
        similarity_matrix = Gene_splice_modeler.compute_similarity_matrix(self.alignments)
        logger.debug("Similarity matrix:\n" + str(similarity_matrix))

        ## build multiple alignment in a hierarchical way
        while len(similarity_matrix) > 1:

            # set diag to -1 to avoid any zero ties w/ self-vals
            for i in range(0,len(alignments)):
                similarity_matrix[ i ][ i ] = -1
            
            ## find best pair
            best_pair_idx = int(numpy.argmax(similarity_matrix))
            num_alignments = len(similarity_matrix)
            best_pair_idx_1 = int(best_pair_idx / num_alignments)
            best_pair_idx_2 = best_pair_idx % num_alignments
            
            ## merge pair into single alignment
            align_a = alignments[ best_pair_idx_1 ]
            align_b = alignments[ best_pair_idx_2 ]

            align_merged = Gene_splice_modeler.merge_alignments(align_a, align_b)
            
            ## recompute matrix
            new_alignment_list = list()
            for i in range(0, len(alignments)):
                if i not in (best_pair_idx_1, best_pair_idx_2):
                    new_alignment_list.append(alignments[ i ])
            new_alignment_list.append(align_merged)

            alignments = new_alignment_list

            logger.debug("\nUpdated alignments:\n" + str(alignments))
            
            similarity_matrix = Gene_splice_modeler.compute_similarity_matrix(alignments)
            logger.debug("Similarity matrix:\n" + str(similarity_matrix))


        if len(alignments) > 1:
            raise RuntimeError("Error, should only have one alignment but have {} alignments after merge".format(len(alignments)))
        
        return alignments[0]


    @staticmethod
    def compute_similarity_matrix(alignments_list):
        """
        similarity matrix indicates number of shared nodes between each pair of isoforms.
        """
        
        num_alignments = len(alignments_list)
        sim_matrix = numpy.zeros( (num_alignments, num_alignments), dtype='int_' )

        for i in range(0, num_alignments-1):
            align_i = alignments_list[i]
            for j in range(i+1, num_alignments):
                align_j = alignments_list[j]

                common_nodes = Node_alignment.Node_alignment.compute_number_common_nodes(align_i, align_j)
                num_common_nodes = len(common_nodes)

                sim_matrix[ i ][ j ] = num_common_nodes
                


        return sim_matrix
        

    @staticmethod
    def merge_alignments(align_a, align_b):
        """
        Computes a mismatch-free multiple alignment (just matches and gaps) between two Node_alignment objects

        returns single Node_alignment object containing the contents of aligned align_a and align_b as aligned.
        
        """
        
        logger.debug("Merging alignments {} and {}".format(align_a, align_b))

        ## ensure the transcripts are disjoint
        transcript_names_align_A = set(align_a.get_transcript_names())
        transcript_names_align_B = set(align_b.get_transcript_names())

        if not set.isdisjoint(transcript_names_align_A, transcript_names_align_B):
            raise RuntimeError("Error, transcripts in alignments to merge are not disjoint: {} and {}".format(transcript_names_align_A, transcript_names_align_B))

        
        width_a = align_a.width()
        width_b = align_b.width()

        # do global alignments w/o penalizing end gaps
        dp_matrix = DP_matrix.DP_matrix.build_DP_matrix(width_a, width_b)

        # put align B across top (cols) and align A at side (row)
        # init the matrix zero rows
        for i in range(1, width_a+1):
            dp_matrix[ i ][ 0 ]['bt'] = 'DEL_B' # UP
        for j in range(1, width_b+1):
            dp_matrix[ 0 ][ j ]['bt'] = 'DEL_A' # LEFT
        
        # score the DP matrix
        for i in range(1, width_a+1):
            for j in range(1, width_b+1):

                score_cell_match = Gene_splice_modeler.get_match_score(align_a, i-1, align_b, j-1) # score matrix is 1-based, align is 0-based
                
                score_diag = dp_matrix[ i-1 ][ j-1 ]['score'] + score_cell_match

                score_del_a = dp_matrix[ i ][ j-1 ]['score']

                score_del_b = dp_matrix[ i-1 ][ j ]['score']


                if score_cell_match > 0 and score_diag >= score_del_a and score_diag >= score_del_b:
                    dp_matrix[ i ][ j ]['score'] = score_diag
                    dp_matrix[ i ][ j ]['bt'] = 'DIAG'
                elif score_del_a >= score_del_b:
                    dp_matrix[ i ][ j ]['score'] = score_del_a
                    dp_matrix[ i ][ j ]['bt'] = 'DEL_A'
                else:
                    dp_matrix[ i ][ j ]['score'] = score_del_b
                    dp_matrix[ i ][ j ]['bt'] = 'DEL_B'


        #logger.debug("DP_matrix:\n" + DP_matrix.toString(dp_matrix))

        """
        # get max score
        max_score = 0
        max_i = -1
        max_j = -1
        for i in range(0,width_a+1):
            score = dp_matrix[ i ][ width_b ]['score']
            if score > max_score:
                max_score = score
                max_i = i
                max_j = width_b
        for j in range(0, width_b+1):
            score = dp_matrix[ width_a ][ j ]['score']
            if score > max_score:
                max_score = score
                max_i = width_a
                max_j = j
        
        logger.info("found max score {} at position: ({},{})".format(max_score, max_i, max_j))
        """

        # keep as global alignment
        max_i = width_a
        max_j = width_b
        
        # backtrack
        i = max_i
        j = max_j
        all_merged_alignment_nodes_list = list()
        while i > 0 or j > 0:
            score_struct = dp_matrix[ i ][ j ]
            
            nodes_align_a = align_a.get_node_LIST_at_column_pos(i-1) # again, remember align has zero-based coords, whereas dp_matrix is 1-based
            nodes_align_b = align_b.get_node_LIST_at_column_pos(j-1)

            align_nodes = list()
                        
            bt_dir = score_struct['bt']

            #logger.debug("backtrack-dir: " + bt_dir)

            if bt_dir == 'DIAG':
                i -= 1
                j -= 1
                align_nodes = nodes_align_a + nodes_align_b
            

            elif bt_dir == 'DEL_B':   # UP
                i -= 1

                align_nodes += nodes_align_a
                for x in range(0,len(nodes_align_b)):
                    align_nodes.append(None)
            
            elif bt_dir == 'DEL_A':  # LEFT
                j -= 1

                for x in range(0,len(nodes_align_a)):
                    align_nodes.append(None)
                align_nodes += nodes_align_b

            else:
                raise RuntimeError("bt: ({},{}), bt_dir not defined".format(i,j))

            all_merged_alignment_nodes_list.append(align_nodes)

        all_merged_alignment_nodes_list.reverse()
        logger.debug("Merged alignment nodes list: " + str(all_merged_alignment_nodes_list) )        

        # prep merged alignment obj
        merged_transcript_name_list = align_a.get_transcript_names() + align_b.get_transcript_names()
        node_obj_matrix = list()
        # interate through each node list, reorganize into a matrix
        for i in range(0,len(merged_transcript_name_list)):
            row = list()
            for node_obj_list in all_merged_alignment_nodes_list:
                row.append(node_obj_list[i])
            node_obj_matrix.append(row)


        logger.debug("merged alignment node matrix:\n" + str(node_obj_matrix))

        merged_alignment_obj = Node_alignment.Node_alignment(align_a.get_gene_id(), merged_transcript_name_list, node_obj_matrix)

        logger.debug("merged alignment obj:\n" + str(merged_alignment_obj))

        #sys.exit(1) # DEBUG
        
        return merged_alignment_obj

    
                
    @staticmethod
    def get_match_score(align_a, idx_a, align_b, idx_b):
        """
        just determines if indices in two transcripts have the same node identifier
        """
        
        node_set_a = align_a.get_node_set_at_column_pos(idx_a)
        node_set_b = align_b.get_node_set_at_column_pos(idx_b)
    
        node_set_a = Node_alignment.Node_alignment.get_node_loc_ids(node_set_a)
        node_set_b = Node_alignment.Node_alignment.get_node_loc_ids(node_set_b)
            
        if (set.intersection(node_set_a, node_set_b)):
            return 1 # match
        else:
            return 0 # no match


    @staticmethod
    def write_malign(gene_name, malign_dict, ofh, align_width=100):
        """
        writes the multiply aligned isoform sequences to an output filehandle
        """
        
        transcript_names = malign_dict.keys()

        alignment_length = len(malign_dict[ transcript_names[ 0 ] ])

        align_start = 0

        align_text = ""

        while align_start < alignment_length:
            for transcript_name in transcript_names:
                align_region = malign_dict[ transcript_name ][ align_start : min(alignment_length, align_start + align_width) ]
                align_text += transcript_name + "\t" + align_region + "\n"
            align_text += "\n" # spacer between alignment blocks
            align_start += align_width

        ofh.write("// {}\n\n{}\n".format(gene_name, align_text))


