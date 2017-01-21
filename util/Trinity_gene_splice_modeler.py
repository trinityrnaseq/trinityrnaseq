#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging
import argparse
import collections

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


class Node:

    node_cache = dict()

    # class method
    def get_node(self, transcript_name, loc_node_id, node_seq):
        node_id = get_node_id(transcript_name, loc_node_id, node_seq)
        if node_id in node_cache:
            node_obj = node_cache[ node_id ]
            if node_obj.seq != node_seq:
                raise RuntimeException("Error, have conflicting node sequences for node_id: {}".format(node_id)
            return node_obj
        else:
            # instantiate a new one
            node_obj = Node(transcript_name, loc_node_id, node_seq)
            node_cache[ node_id ] = node_obj
            return node_obj

    # class method
    def get_node_id(self, transcript_name, loc_node_id, node_seq):
        node_id = "::".join(transcript_name, loc_node_id, len(node_seq))
        return node_id


    def __init__(self, transcript_name, loc_node_id, node_seq):
        self.transcript_name = transcript_name
        self.loc_node_id = loc_node_id
        self.seq = node_seq
        self.len = len(node_seq)

        self.prev = set()
        self.next = set()

    def add_next_node(self, next_node_obj):
        self.next.add(next_node_obj)

    def add_prev_node(self, prev_node_obj):
        self.prev.add(prev_node_obj)



class Node_path:

    transcript_name = ""
    node_obj_list = list()

    def __init__(self, transcript_name, path_string, sequence):
        self.transcript_name = transcript_name
        int_node_list = re.findall("\d+", path_string)

        obj_node_list = list()
        for int_value in int_node_list:
            int_value = int(int_value)
            node_obj = node(transcript_name, int_value)

            obj_node_list.append(node_obj)


class Trinity_fasta_parser:

    trinity_gene_to_isoform_seqs = collections.defaultdict(list)

    def __init__(self, trinity_fasta_filename):

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
            raise RuntimeException("Error, cannot parse accession from header: {}".format(header))
        accession = m.group(1)

        m = re.search("path=\[([^\]]+)\]", header)
        if not m:
            raise RuntimeExcpetion("Error, cannot parse path info from header of line: {}".format(header))
        
        path_str = m.group(1)
        
        # get gene ID
        gene_id = re.sub("_i\d+$", "", accession)
        if gene_id == accession:
            raise RuntimeException("Error, couldn't remove isoform ID from Trinity accession: {}".format(accession))

        isoform_list = trinity_gene_to_isoform_seqs[ gene_id ]

        iso_struct = { 'transcript_name' : accession,
                       'path' : path_str,
                       'seq' : sequence }

        isoform_list.append(iso_struct)
        
    


def main():

    parser = argparse.ArgumentParser(description="Converts Trinity Isoform structures into a single gene structure representation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--trinity_fasta", dest="trinity_fasta", type=str, default="", required=True, help="Trinity.fasta file")

    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)      


    gene_to_isoform_paths = parse_gene_to_iso_paths(args.trinity_fasta)



 
####################
 
if __name__ == "__main__":
    main()
