#!/usr/bin/env python

import sys, os, re
import subprocess

usage = "\n\n\tusage: " + sys.argv[0] + " pct_align.stats max_pct_aligned transcript.fasta\n\n"

if len(sys.argv) < 4:
    print >>sys.stderr, usage
    sys.exit(1)

pct_align_stats_filename = sys.argv[1]
max_pct_aligned = float(sys.argv[2])
transcripts_fasta_filename = sys.argv[3]


def main():

    if not os.path.exists(transcripts_fasta_filename + ".fai"):
        if os.system("samtools faidx " + transcripts_fasta_filename):
            raise Exception("Error, cannot make fai file") 

    for line in open(pct_align_stats_filename):
        #print line
        line = line.rstrip()
        x = line.split("\t")
        trans_acc = x[0]
        pct_aligned = float(x[3])

        if pct_aligned < max_pct_aligned:
            
            cmd = "samtools faidx " + transcripts_fasta_filename + " \'" + trans_acc + "\'"
            fasta_seq = subprocess.check_output(cmd, shell=True)
            sys.stdout.write(fasta_seq)



    sys.exit(0)



main()
