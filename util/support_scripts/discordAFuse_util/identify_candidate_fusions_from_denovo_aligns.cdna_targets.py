#!/usr/bin/env python

import sys, os, re

usage = "\n\n\tusage: " + sys.argv[0] + " chimera.trans.align.gff3 cdna.fasta\n\n"

if len(sys.argv) < 3:
    print >>sys.stderr, usage
    sys.exit(1)

chimera_trans_align_gff3 = sys.argv[1]
cdna_fasta_filename = sys.argv[2]



def main():
    
    headers_file = cdna_fasta_filename + ".headers"
    if not os.path.exists(headers_file):
        print >>sys.stderr, "Error, cannot locate headers file: " + headers_file
        sys.exit(2)

    trans_to_gene_hash = parse_trans_to_gene_from_header(headers_file)

    trans_to_gene_pairs = parse_trans_to_gene_pairs(chimera_trans_align_gff3, trans_to_gene_hash)

    for trans_acc in trans_to_gene_pairs:
        gene_hits_array = trans_to_gene_pairs[trans_acc]
        #if gene_hits_array[0][0] == None or gene_hits_array[1][0] == None: continue

        gene_1_name = gene_hits_array[0][0]
        gene_1_orient = gene_hits_array[0][1]

        gene_2_name = gene_hits_array[1][0]
        gene_2_orient = gene_hits_array[1][1]

        #print "\t".join(str(x) for x in [trans_acc, gene_1_name, gene_1_orient, gene_2_name, gene_2_orient])

        #continue
        
        if gene_1_orient != gene_2_orient: continue

        if (gene_1_orient, gene_2_orient) == ('-','-'):
            # swap them
            (gene_1_name, gene_2_name) = (gene_2_name, gene_1_name)

        #print "\t".join([trans_acc, gene_1_name, gene_2_name])
        print gene_1_name + "--" + gene_2_name + "\t" + trans_acc

    

    sys.exit(0)



def parse_trans_to_gene_pairs(chimera_gff3, trans_to_gene_hash):

    trans_to_gene_pairs = {}

    for line in open(chimera_gff3):
        if re.match("^\#", line): continue
        if not re.match("\w", line): continue
        #print "Line: " + line

        x = line.split("\t")

        info = x[8]
        trans_acc = x[0]
        orient = x[6]
        
        path_num = re.search(".path(\d);", info).group(1)
        path_num = int(path_num)
        path_num -= 1

        denovo_trans_id = re.search("Target=(\S+)", info).group(1)

        gene_name = trans_to_gene_hash[trans_acc]

        if denovo_trans_id not in trans_to_gene_pairs:
            trans_to_gene_pairs[denovo_trans_id] = [ [None,None], [None,None] ]
            
        trans_to_gene_pairs[denovo_trans_id][path_num] = [gene_name, orient]
    


    return trans_to_gene_pairs



def parse_trans_to_gene_from_header(headers_file):

    trans_to_gene_hash = {}

    for line in open(headers_file):
        line = line.rstrip()
        x = re.split("\s+", line)
        trans_id = x[0]
        gene_symbol = x[2] if len(x) > 2 else x[1]

        trans_to_gene_hash[trans_id] = gene_symbol

    return trans_to_gene_hash




####
main()
