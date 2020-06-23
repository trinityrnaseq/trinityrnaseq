#!/usr/bin/env python3
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

# add local py lib
sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "pylib"]))

import TGraph
import Trinity_fasta_parser
import Gene_splice_modeler
import Splice_model_refiner
import Node_path
import TGLOBALS

def main():

    parser = argparse.ArgumentParser(
        description="Converts Trinity Isoform structures into a single gene structure representation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--trinity_fasta", dest="trinity_fasta", type=str,
                        default="", required=True, help="Trinity.fasta file")

    parser.add_argument("--out_prefix", dest="out_prefix", type=str,
                    default="trinity_genes", required=False, help="output prefix for fasta and gtf outputs")
    parser.add_argument("--incl_malign", dest="malign", action="store_true",
                        default=False, help="include multiple alignment formatted output file")
    
    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")
    parser.add_argument("--verbose", required=False, action="store_true", default=False, help="verbose mode")

    parser.add_argument("--no_squeeze", required=False, action="store_true",
                        default=False, help="don't merge unbranched stretches of node identifiers")

    parser.add_argument("--no_refinement", required=False, action="store_true",
                        default=False, help="don't refine splice graph by further collapsing allelic variants")
    
    parser.add_argument("--incl_cdna", required=False, action="store_true",
                        default=False, help="rewrite Trinity fasta file using simplified graph structure")

    parser.add_argument("--incl_dot", required=False, action="store_true",
                        default=False, help="include dot file for gene graph (*warning* single dot file per gene!! use sparingly)")

    parser.add_argument("--restrict_gene_id", required=False, dest='restrict_gene_id', type=str,
                        default=None, help="only process this gene")
    
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        TGLOBALS.DEBUG = True
        args.verbose = True

    logger.info("-parsing Trinity fasta file: {}".format(args.trinity_fasta))
    trin_parser = Trinity_fasta_parser.Trinity_fasta_parser(args.trinity_fasta)

    logger.info("-organizing gene/isoform data")
    gene_to_isoform_info = trin_parser.get_trinity_gene_to_isoform_info()

    out_fasta_filename = args.out_prefix + ".fasta"
    out_gtf_filename = args.out_prefix + ".gtf"
    out_trinity_fa_filename = args.out_prefix + ".transcripts.fa"
    out_malign_filename = args.out_prefix + ".malign"
    

    ofh_fasta = open(out_fasta_filename, 'w')
    ofh_gtf = open(out_gtf_filename, 'w')
    ofh_trinity_fa = None
    ofh_malign = None

    if args.malign:
        ofh_malign = open(out_malign_filename, 'w')

    if args.incl_cdna:
        ofh_trinity_fa = open(out_trinity_fa_filename, 'w')



    supertranscript_start_time = time.time()
    
    ## examine the alt spliced isoforms.

    logger.info("-computing supertranscripts")
    num_genes = len(gene_to_isoform_info.keys())
    gene_counter = 0
    
    for gene_name in gene_to_isoform_info:

        if args.restrict_gene_id and args.restrict_gene_id != gene_name:
            continue

        iso_struct_list = gene_to_isoform_info[ gene_name ]

        gene_counter += 1

        tgraph = TGraph.TGraph(gene_name)

        # convert to Node_path objects
        node_path_obj_list = list()
        for iso_struct in iso_struct_list:
            n_path = Node_path.Node_path(tgraph, iso_struct['transcript_name'], iso_struct['path'], iso_struct['seq'])
            node_path_obj_list.append(n_path)
            #print(str(n_path))

        # adjust for unique prefix yet shared suffix of FST nodes
        node_path_obj_list = Node_path.Node_path.adjust_for_fst_nodes(tgraph, node_path_obj_list)

        
        # generate multiple alignment

        start_time = time.time()

        logger.debug("\n\n##################\n############## PROCESSING GENE: {}\n###########################\n\n".format(gene_name))
        logger.info("# Processing Gene: {} having {} isoforms".format(gene_name, len(node_path_obj_list)))
        
        gene_splice_modeler = Gene_splice_modeler.Gene_splice_modeler(gene_name, node_path_obj_list)
        
        splice_model_alignment = gene_splice_modeler.build_splice_model()


        if args.verbose:
            logger.info("Final splice_model_alignment for Gene {} :\n{}\n".format(gene_name, str(splice_model_alignment)))


        ################################
        ## Splice graph refinement stage

        logger.info("Splice graph refinement underway")
        squeezed_splice_model = splice_model_alignment

        
        if args.no_squeeze:
            logger.info("--no_squeeze set, so not squeezing structure")
            #squeezed_splice_model = Splice_model_refiner.refine_alignment(squeezed_splice_model, reset_node_ids=True) # True important here... can have repeat nodes
        else:
            # SQUEEZE
            squeezed_splice_model = splice_model_alignment.squeeze()
            logger.debug("Squeezed splice model for Gene {}:\n{}\n".format(gene_name, str(squeezed_splice_model)))

        if not args.no_refinement:
            squeezed_splice_model = Splice_model_refiner.refine_alignment(squeezed_splice_model,
                                                                          reset_node_ids=True, # True important here... can have repeat nodes
                                                                          max_burr_length=25,
                                                                          max_bubble_pop_length=25)
            
            logger.debug("After seq-refinement of splice model for Gene {}:\n{}\n".format(gene_name, str(squeezed_splice_model)))
            squeezed_splice_model = Splice_model_refiner.remove_redundant_paths(squeezed_splice_model)
            logger.debug("After removing redundant paths of splice model for Gene {}:\n{}\n".format(gene_name, str(squeezed_splice_model)))
            
            if not args.no_squeeze:
                # re-squeeze
                squeezed_splice_model = squeezed_splice_model.squeeze()
                logger.debug("Post-refinement and final squeeze of splice model for Gene {}:\n{}\n".format(gene_name, str(squeezed_splice_model)))


        # remove redundant sequence entries (created during the refinement process)
        if squeezed_splice_model.remove_redundant_sequences():
            # resqueeze
            squeezed_splice_model = squeezed_splice_model.squeeze()
        
                
        # done with splice graph refinement.
                
        squeezed_splice_model.reassign_node_loc_ids_by_align_order()
        
        
        (gene_seq, gtf_txt, trinity_fa_text, malign_dict) = squeezed_splice_model.to_gene_fasta_and_gtf(gene_name)

        if args.incl_dot:
            squeezed_splice_model.to_splice_graph(gene_name).draw_graph("{}.dot".format(gene_name))
            write_mfa(malign_dict, "{}.mfa".format(gene_name))
        
        ofh_fasta.write(">{}\n{}\n".format(gene_name, gene_seq))
        ofh_gtf.write(gtf_txt + "\n")

        if args.incl_cdna:
            ofh_trinity_fa.write(trinity_fa_text)

        if args.malign and len(node_path_obj_list) > 1:
            Gene_splice_modeler.Gene_splice_modeler.write_malign(gene_name, malign_dict, ofh_malign)


        runtime = time.time() - start_time
        if runtime > 0.1 or args.debug:
            pct_done = float(gene_counter)/num_genes * 100
            logger.info("Exec Time for Gene {}: {:.3f} s, total pct done: {:.2f}%\n".format(gene_name, runtime, pct_done))
        

        if args.restrict_gene_id and args.restrict_gene_id == gene_name:
            break
    

    ofh_fasta.close()
    ofh_gtf.close()
    if args.incl_cdna:
        ofh_trinity_fa.close()
    if args.malign:
        ofh_malign.close()

    supertranscript_end_time = time.time()
    runtime_minutes = (supertranscript_end_time - supertranscript_start_time) / 60.0

    logger.info("Done.  Total runtime: {:.1f} min\n\n".format(runtime_minutes))
    

    sys.exit(0)


def write_mfa(malign_dict, filename):

    ofh = open(filename, 'w')
    for acc in malign_dict:
        aligned_seq = malign_dict[acc]
        ofh.write(">{}\n{}\n".format(acc, aligned_seq))

    ofh.close()


 
####################
 
if __name__ == "__main__":
    main()
