#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import logging
import argparse
import gzip
from collections import defaultdict

def main():

    parser = argparse.ArgumentParser(description="validate fastqs", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    
    parser.add_argument("--left_fq", required=True, type=str, help="left fastq file")

    parser.add_argument("--right_fq", required=False, type=str, help="right fastq file")
    
    args = parser.parse_args()

    left_fq = args.left_fq
    right_fq = args.right_fq




    left_fq_iterator = fastq_iterator(left_fq)
    right_fq_iterator = None
    if right_fq:
        right_fq_iterator = fastq_iterator(right_fq)



    counter = 0
    
    for left_read_tuple in left_fq_iterator:
        left_readname, left_readseq, left_L3, left_quals = left_read_tuple
        left_readseq_len = len(left_readseq)
        left_quals_len = len(left_quals)

        counter += 1
        if counter % 10000 == 0:
            sys.stderr.write(f"\r[{counter}]  ")

        
        assert left_readseq_len == left_quals_len, f"Error, left seqlen and quals len dont match:\n{left_readname}\n{left_readseq}\n{left_L3}\n{left_quals}\n"
        
        
        if right_fq_iterator is not None:

            right_read_tuple = next(right_fq_iterator)
            right_readname, right_readseq, right_L3, right_quals = right_read_tuple
            right_readseq_len = len(right_readseq)
            right_quals_len = len(right_quals)
            

            assert right_readseq_len == right_quals_len, f"Error, right seqlen and quals len dont match:\n{right_readname}\n{right_readseq}\n{right_L3}\n{right_quals}\n"

            assert core_readname(left_readname) == core_readname(right_readname), f"Error, fastq pair read names aren't consistent: {left_readname} vs. {right_readname}"
            

    print("fastq(s) validate")
    
    sys.exit(0)



def core_readname(readname):

    core_readname = readname.split(" ")[0]

    core_readname = re.sub("/[12]$", "", core_readname)

    return core_readname



def fastq_iterator(fastq_filename):

    if re.search(".gz$", fastq_filename):
        fh = gzip.open(fastq_filename, 'rt', encoding='utf-8')
    else:
        fh = open(fastq_filename, 'rt', encoding='utf-8')

    have_records = True
    while have_records:
        readname = next(fh).rstrip()
        readseq = next(fh).rstrip()
        L3 = next(fh).rstrip()
        quals = next(fh).rstrip()

        yield (readname, readseq, L3, quals)

        if not readname:
            break

    return
        



if __name__=='__main__':
    main()
