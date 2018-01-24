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


class Trinity_fasta_parser:
    """
    Parses a Trinity.fasta file and stores the transcript name, sequence, and node path info.

    Instance member:

        trinity_gene_to_isoform_seqs : (defaultdict(list)) stores key,val of transcript_name,path_struct

        where path_struct has structure:
             {
                 'transcript_name' : accession,
                 'path' : path_str,
                 'seq' : sequence
             }
    """
    
    def __init__(self, trinity_fasta_filename):

        self.trinity_gene_to_isoform_seqs = collections.defaultdict(list) #instantiate member 

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
            raise RuntimeError("Error, cannot parse path info from header of line: {}".format(header))
        
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


