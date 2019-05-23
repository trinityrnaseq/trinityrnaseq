#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging
import argparse
import collections

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__file__)
 
def main():

    parser = argparse.ArgumentParser(description="estimates gene length as isoform lengths weighted by TPM expression values")
    
    parser.add_argument("--gene_trans_map", dest="gene_trans_map_file", type=str, default="",
                        required=True, help="gene-to-transcript mapping file, format: gene_id(tab)transcript_id")

    parser.add_argument("--trans_lengths", dest="trans_lengths_file", type=str, required=True,
                        help="transcript length file, format: trans_id(tab)length")

    parser.add_argument("--TPM_matrix", dest="TPM_matrix_file", type=str, default="",
                        required=True, help="TPM expression matrix")

    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")


    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
    

    trans_to_gene_id_dict = parse_gene_trans_map(args.gene_trans_map_file)

    trans_lengths_dict = parse_trans_lengths_file(args.trans_lengths_file)

    trans_to_TPM_vals_dict = parse_TPM_matrix(args.TPM_matrix_file)

    weighted_gene_lengths = compute_weighted_gene_lengths(trans_to_gene_id_dict,
                                                          trans_lengths_dict,
                                                          trans_to_TPM_vals_dict)

    print("#gene_id\tlength")
    for gene_id,length in weighted_gene_lengths.items():
        print("\t".join([gene_id,str(length)]))


    sys.exit(0)
    


def compute_weighted_gene_lengths(trans_to_gene_id_dict, trans_lengths_dict, trans_to_TPM_vals_dict):

    gene_id_to_trans_list = collections.defaultdict(list)

    gene_id_to_length = {}

    pseudocount = 1

    for trans_id,gene_id in trans_to_gene_id_dict.items():
        gene_id_to_trans_list[gene_id].append(trans_id)

    for gene_id,trans_list in gene_id_to_trans_list.items():

        if len(trans_list) == 1:

            gene_id_to_length[gene_id] = trans_lengths_dict[ trans_list[0] ]
        else:
            
            sum_length_x_expr = 0
            sum_expr = 0

            trans_expr_lengths = []

            for trans_id in trans_list:
                trans_len = trans_lengths_dict[trans_id]
                expr_vals = trans_to_TPM_vals_dict[trans_id]
                trans_sum_expr = sum(expr_vals) + pseudocount

                trans_expr_lengths.append((trans_len, trans_sum_expr))

                sum_length_x_expr += trans_sum_expr * trans_len
                sum_expr += trans_sum_expr


            weighted_gene_length = sum_length_x_expr / sum_expr
            gene_id_to_length[gene_id] = int(round(weighted_gene_length))

            logger.debug("Computing weighted length of {0}: {1} => {2}".format(gene_id,
                                                                            trans_expr_lengths,
                                                                            weighted_gene_length))
            
    return gene_id_to_length
            


def parse_TPM_matrix(TPM_matrix_file):

    trans_to_TPM_vals_dict = {}
    
    with open(TPM_matrix_file) as f:
        header = next(f)
        for line in f:
            line = line.rstrip()
            vals = line.split("\t")
            trans_id = vals[0]
            expr_vals_list = vals[1:]
            expr_vals_list = [float(x) for x in expr_vals_list]
            trans_to_TPM_vals_dict[trans_id] = expr_vals_list

    return trans_to_TPM_vals_dict


def parse_trans_lengths_file(trans_lengths_file):

    trans_id_to_length = {}

    with open(trans_lengths_file) as f:
        for line in f:
            line = line.rstrip()
            if line[0] == '#':
                continue
            
            (trans_id, length) = line.split("\t")

            if re.match("^\d+$", length):
                trans_id_to_length[trans_id] = int(length)
            else:
                print("Warning - ignoring line: [{0}] since not parsing length value as number".format(line), file=sys.stderr)

    return trans_id_to_length



def parse_gene_trans_map(gene_trans_map_file):

    trans_to_gene_id = {}

    with open(gene_trans_map_file) as f:
        for line in f:
            line = line.rstrip()
            (gene_id, trans_id) = line.split("\t")

            trans_to_gene_id[trans_id] = gene_id;
        

    return trans_to_gene_id


 
####################
 
if __name__ == "__main__":
    main()
