#!/usr/bin/env python3
import sys
import re
import numpy as np

#regex for cigar
re_split_cigar  = re.compile('(\d+[MDI])(?=[0-9]|$)')

#SAM INDEXES
idx_qname = 0
idx_flag  = 1
idx_rname = 2
idx_pos   = 3
idx_mapq  = 4
idx_cigar = 5
idx_rnext = 6
idx_pnext = 7
idx_tlen  = 8
idx_seq   = 9
idx_qual  = 10

def write_err(msg, exit=False, status=1):
    sys.stderr.write(msg)
    if exit:
        sys.exit(status)

if __name__ == '__main__':
    chrlen   = int(sys.argv[1])
    _chr     = sys.argv[2]
    infile   = sys.argv[3]
    outfile  = sys.argv[4]
    nosingle = int(sys.argv[5])
    minins   = int(sys.argv[6])
    maxins   = int(sys.argv[7]) 

    #open input sam
    if infile == "-":
        insam = sys.stdin

    else:
        try:
            insam = open(infile, 'r')
        except IOError:
            write_err("ERROR: Unable to read input sam file\n", exit=True)

    #get file for writing
    if outfile == "-":
        outwig = sys.stdout
        
    else:
        try:
            outwig = open(outfile, 'w')
        except IOError:
            write_err("ERROR: Unable to open output file\n", exit=True)

    #start generating wig
    wig = np.zeros(chrlen, dtype = np.uint64)

    for row in insam:
        entry = row.split('\t')
        start = int(entry[idx_pos]) - 1
        tlen  = int(entry[idx_tlen])
        atlen = abs(tlen)

        #entry = gen_sam(row)
        #start = entry.pos - 1
        #atlen = abs(entry.tlen)
        if atlen > 0 and atlen >= minins and atlen <= maxins:
            if tlen > 0:
                wig[start:(start + tlen)] += 1
            else:
                continue
        else:
            cigars = re.split(re_split_cigar, entry[idx_cigar])
            cigars = [x for x in cigars if x != '']
            _len   = 0

            for c in cigars:
                if 'M' in c:
                    _len += int(c.split('M')[0])
                elif 'D' in c:
                    _len += int(c.split('D')[0])
                elif 'I' in c:
                    continue
                    _len += int(c.split('I')[0])
                else:
                    continue
            wig[start:(start + _len)] += 1


    outwig.write("variableStep chrom=" + _chr + "\n")
    for i in range(0, chrlen):
        outwig.write(str(i + 1) + "\t" + str(wig[i]) + "\n")


