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
        self.transcript_name = transcript_name
        self.loc_node_id = loc_node_id
        self.seq = node_seq
        self.len = len(node_seq)

        self.prev = set()
        self.next = set()


    @classmethod
    def get_node(self, transcript_name, loc_node_id, node_seq):
        node_id = get_node_id(transcript_name, loc_node_id, node_seq)
        if node_id in node_cache:
            node_obj = node_cache[ node_id ]
            if node_obj.seq != node_seq:
                raise RuntimeError("Error, have conflicting node sequences for node_id: {}".format(node_id))
            return node_obj
        else:
            # instantiate a new one
            node_obj = Node(transcript_name, loc_node_id, node_seq)
            node_cache[ node_id ] = node_obj
            return node_obj

    @staticmethod
    def get_node_id(transcript_name, loc_node_id, node_seq):
        node_id = "::".join([transcript_name, loc_node_id, str(len(node_seq))])
        return node_id

    
    def get_loc_id(self):
        return self.loc_node_id
    
    def get_seq(self):
        return self.seq

    def get_transcript_name(self):
        return self.transcript_name
        

    def add_next_node(self, next_node_obj):
        self.next.add(next_node_obj)

    def add_prev_node(self, prev_node_obj):
        self.prev.add(prev_node_obj)

    def __repr__(self):
        #return(self.get_node_id(self.transcript_name, self.loc_node_id, self.seq))
        return(self.loc_node_id)
        
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



class Node_path:

    def __init__(self, transcript_name, path_string, sequence):
        self.transcript_name = transcript_name
        self.node_obj_list = list()

        node_descr_list = re.findall("\d+:\d+\-\d+", path_string)

        obj_node_list = list()
        for node_descr in node_descr_list:
            (loc_node_id, node_coord_range) = node_descr.split(":")
            (lend,rend) = node_coord_range.split("-")
            lend = int(lend)
            rend = int(rend)
            
            node_obj = Node(transcript_name, loc_node_id, sequence[lend:rend+1]) # coords in path were already zero-based

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
        

class Trinity_fasta_parser:

    def __init__(self, trinity_fasta_filename):

        self.trinity_gene_to_isoform_seqs = collections.defaultdict(list)

        with open(trinity_fasta_filename) as fh:
            header = ""
            sequence = ""
            for line in fh:
                line = line.rstrip()
                if line[0] == '>':
                    # got header line
                    # process earlier entry
                    if header != "" and sequence != "":
                        self.add_trinity_seq_entry(header, sequence)
                    # init new entry                        
                    header = line
                    sequence = ""
                else:
                    # sequence line
                    sequence += line
            # get last one
            if sequence != "":
                self.add_trinity_seq_entry(header, sequence)


    def add_trinity_seq_entry(self, header, sequence):
        """
        entry looks like so:
        >TRINITY_DN16_c0_g1_i2 len=266 path=[1:0-48 27:49-49 28:50-50 27:51-51 28:52-52 27:53-53 28:54-54 27:55-55 28:56-56 27:57-57 28:58-58 27:59-59 28:60-60 27:61-61 29:62-265] [-1, 1, 27, 28, 27, 28, 27, 28, 27, 28, 27, 28, 27, 28, 27, 29, -2]
        CTGTTGTGTGGGGGGTGCGCTTGTTTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC
        TCAAGTTGATTCCTCCATGTTGCTTTACAGAGACCTGCCAACTACCCAGGAATGTAAAAG
        CATTCATAGTATTTGTCTAGTAGAGATGCTGTATGAAAAATGCCAAAACCAAAAAGAGAA
        AGAAGGAAAGAGAGATAGATAGATGACATAGATGACGGATGGATGGGTGGGTGGGTGGAT
        GGATGGATGGATGGATGGAGGGGGGC
        """

        m = re.search("^>(\S+)", header)
        if not m:
            raise RuntimeError("Error, cannot parse accession from header: {}".format(header))
        accession = m.group(1)

        m = re.search("path=\[([^\]]+)\]", header)
        if not m:
            raise RuntimeExcpetion("Error, cannot parse path info from header of line: {}".format(header))
        
        path_str = m.group(1)
        
        # get gene ID
        gene_id = re.sub("_i\d+$", "", accession)
        if gene_id == accession:
            raise RuntimeError("Error, couldn't remove isoform ID from Trinity accession: {}".format(accession))

        isoform_list = self.trinity_gene_to_isoform_seqs[ gene_id ]

        iso_struct = { 'transcript_name' : accession,
                       'path' : path_str,
                       'seq' : sequence }

        isoform_list.append(iso_struct)
        
    def get_trinity_gene_to_isoform_info(self):
        return self.trinity_gene_to_isoform_seqs


class Node_alignment:

    GAP = None

    def __init__(self, transcript_name_list, node_obj_matrix):
        self.transcript_names = transcript_name_list
        self.aligned_nodes = node_obj_matrix

    def get_transcript_names(self):
        return self.transcript_names

    def get_aligned_nodes(self):
        return self.aligned_nodes

    @staticmethod
    def get_single_seq_node_alignment(transcript_name, path_obj):
        node_list = list()
        for node_obj in  path_obj.get_path():
            node_list.append(node_obj)

        self = Node_alignment([transcript_name], [node_list])

        return self


    @staticmethod
    def compute_number_common_nodes(align_A, align_B):
        node_set_a = Node_alignment.get_node_set(align_A)
        node_set_b = Node_alignment.get_node_set(align_B)

        node_set_a = Node_alignment.get_node_loc_ids(node_set_a)
        node_set_b = Node_alignment.get_node_loc_ids(node_set_b)
                
        common_nodes = set.intersection(node_set_a, node_set_b)

        return common_nodes


    @staticmethod
    def get_node_loc_ids(node_set):
        loc_ids_set = set()
        for node in node_set:
            loc_id = node.get_loc_id()
            loc_ids_set.add(loc_id)

        return loc_ids_set
    

    @staticmethod
    def get_node_set(align_obj):
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
                    node_to_add = Node.merge_nodes(blocked_nodes)
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
        

class Gene_splice_modeler:

    def __init__(self, node_path_obj_list):

        ## create alignments
        self.alignments = list()

        logger.debug("Gene_splice_modeler inputs: {}".format(node_path_obj_list))
        
        for node_path_obj in node_path_obj_list:
            transcript_name = node_path_obj.get_transcript_name()
            alignment_obj = Node_alignment.get_single_seq_node_alignment(transcript_name, node_path_obj)

            self.alignments.append(alignment_obj)

            #print(alignment_obj)
            
    def build_splice_model(self):

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
        num_alignments = len(alignments_list)
        sim_matrix = numpy.zeros( (num_alignments, num_alignments), dtype='int_' )

        for i in range(0, num_alignments-1):
            align_i = alignments_list[i]
            for j in range(i+1, num_alignments):
                align_j = alignments_list[j]

                common_nodes = Node_alignment.compute_number_common_nodes(align_i, align_j)
                num_common_nodes = len(common_nodes)

                sim_matrix[ i ][ j ] = num_common_nodes
                


        return sim_matrix
        

    @staticmethod
    def merge_alignments(align_a, align_b):

        logger.debug("Merging alignments {} and {}".format(align_a, align_b))

        ## ensure the transcripts are disjoint
        transcript_names_align_A = set(align_a.get_transcript_names())
        transcript_names_align_B = set(align_b.get_transcript_names())

        if not set.isdisjoint(transcript_names_align_A, transcript_names_align_B):
            raise RuntimeError("Error, transcripts in alignments to merge are not disjoint: {} and {}".format(transcript_names_align_A, transcript_names_align_B))

        
        width_a = align_a.width()
        width_b = align_b.width()

        # do global alignments w/o penalizing end gaps
        dp_matrix = DP_matrix.build_DP_matrix(width_a, width_b)

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

        merged_alignment_obj = Node_alignment(merged_transcript_name_list, node_obj_matrix)

        logger.debug("merged alignment obj:\n" + str(merged_alignment_obj))

        #sys.exit(1) # DEBUG
        
        return merged_alignment_obj

    
                
    @staticmethod
    def get_match_score(align_a, idx_a, align_b, idx_b):
        
        node_set_a = align_a.get_node_set_at_column_pos(idx_a)
        node_set_b = align_b.get_node_set_at_column_pos(idx_b)
    
        node_set_a = Node_alignment.get_node_loc_ids(node_set_a)
        node_set_b = Node_alignment.get_node_loc_ids(node_set_b)
            
        if (set.intersection(node_set_a, node_set_b)):
            return 1 # match
        else:
            return 0 # no match


    @staticmethod
    def write_malign(gene_name, malign_dict, ofh, align_width=100):

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

                
class DP_matrix:

    @staticmethod
    def build_DP_matrix(num_rows, num_cols):
        dp_matrix = list()
        for i in range(0, num_rows+1):
            row = []
            for j in range(0, num_cols+1):
                struct = { 'i' : i,
                           'j' : j,
                           'bt' : -1,
                           'score' : 0,
                           }
                row.append(struct)
            dp_matrix.append(row)
        
        return dp_matrix    


    @staticmethod
    def toString(dp_matrix):
        nrow = len(dp_matrix)
        ncol = len(dp_matrix[0])

        logger.debug("Matrix is {} X {}".format(nrow, ncol))

        ret_text = ""
        for row in dp_matrix:
            str_row = [ str(x['score']) for x in row ]
            ret_text += "\t".join(str_row) + "\n"

        return ret_text


def main():

    parser = argparse.ArgumentParser(description="Converts Trinity Isoform structures into a single gene structure representation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--trinity_fasta", dest="trinity_fasta", type=str, default="", required=True, help="Trinity.fasta file")

    parser.add_argument("--out_prefix", dest="out_prefix", type=str, default="trinity_genes", required=False, help="output prefix for fasta and gtf outputs")
    parser.add_argument("--incl_malign", dest="malign", action="store_true", default=False, help="include multiple alignment formatted output file")
    
    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)      


    trin_parser = Trinity_fasta_parser(args.trinity_fasta)

    gene_to_isoform_info = trin_parser.get_trinity_gene_to_isoform_info()

    out_fasta_filename = args.out_prefix + ".fasta"
    out_gtf_filename = args.out_prefix + ".gtf"
    out_malign_filename = args.out_prefix + ".malign"

    ofh_fasta = open(out_fasta_filename, 'w')
    ofh_gtf = open(out_gtf_filename, 'w')
    ofh_malign = None

    if args.malign:
        open(out_malign_filename, 'w')


    supertranscript_start_time = time.time()
    
    ## examine the alt spliced isoforms.

    num_genes = len(gene_to_isoform_info.keys())
    gene_counter = 0
    
    for gene_name in gene_to_isoform_info:
        iso_struct_list = gene_to_isoform_info[ gene_name ]

        gene_counter += 1
        
        # convert to Node_path objects
        node_path_obj_list = list()
        for iso_struct in iso_struct_list:
            n_path = Node_path(iso_struct['transcript_name'], iso_struct['path'], iso_struct['seq'])
            node_path_obj_list.append(n_path)
            #print(str(n_path))
        
        # generate multiple alignment

        start_time = time.time()
        
        logger.info("Processing Gene: {} having {} isoforms".format(gene_name, len(node_path_obj_list)))

        gene_splice_modeler = Gene_splice_modeler(node_path_obj_list)

        splice_model_alignment = gene_splice_modeler.build_splice_model()

        logger.debug("Final splice_model_alignment for Gene {} :\n{}\n".format(gene_name, str(splice_model_alignment)))

        squeezed_splice_model = splice_model_alignment.squeeze()
        
        logger.debug("Squeezed splice model for Gene {}:\n{}\n".format(gene_name, str(squeezed_splice_model)))
        
        (gene_seq, gtf_txt, malign_dict) = squeezed_splice_model.to_gene_fasta_and_gtf(gene_name)

        ofh_fasta.write(">{}\n{}\n".format(gene_name, gene_seq))
        ofh_gtf.write(gtf_txt + "\n")

        if args.malign and len(node_path_obj_list) > 1:
            Gene_splice_modeler.write_malign(gene_name, malign_dict, ofh_malign)


        runtime = time.time() - start_time
        if runtime > 0.1 or args.debug:
            pct_done = float(gene_counter)/num_genes * 100
            logger.info("Exec Time for Gene {}: {:.3f} s, total pct done: {:.2f}%\n".format(gene_name, runtime, pct_done))
        

    ofh_fasta.close()
    ofh_gtf.close()
    if args.malign:
        ofh_malign.close()

    supertranscript_end_time = time.time()
    runtime_minutes = (supertranscript_end_time - supertranscript_start_time) / 60.0

    logger.info("Done.  Total runtime: {:.1f} min\n\n".format(runtime_minutes))
    

    sys.exit(0)

 
####################
 
if __name__ == "__main__":
    main()
