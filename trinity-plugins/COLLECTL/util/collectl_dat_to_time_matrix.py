#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging
import argparse
import collections
import datetime
import time

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(description="generate time matrix from collectl dat file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--dat", dest="dat_file", type=str, default="", required=True,
                        help="collectl dat file")

    parser.add_argument("--out_prefix", dest="out_prefix", type=str, default="collectl",
                        help="prefix for output files")

    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)      

    
    cpu_usage_matrix = collections.defaultdict(dict)
    memory_usage_matrix = collections.defaultdict(dict)
    times = list()
    prognames = list()

    parse_collectl_dat(args.dat_file, cpu_usage_matrix, memory_usage_matrix, times, prognames)

    output_prefix = args.out_prefix
    cpu_usage_output_filename = output_prefix + ".cpu_usage.matrix"
    mem_usage_output_filename = output_prefix + ".mem_usage.matrix"

    cpu_usage_ofh = open(cpu_usage_output_filename, 'w')
    mem_usage_ofh = open(mem_usage_output_filename, 'w')

    # print headers
    header = "\t" + "\t".join(prognames)
    cpu_usage_ofh.write(header + "\n")
    mem_usage_ofh.write(header + "\n")
    
    for timeval in times:
        memvals = list()
        cpuvals = list()
        
        timeval = str(timeval)

        memvals.append(timeval)
        cpuvals.append(timeval)
        
        for progname in prognames:
            memory_hash = memory_usage_matrix[timeval]
            memory_val = "0"
            if progname in memory_hash:
                memory_val = str(memory_hash[progname])
            memvals.append(memory_val)

            cpu_val = "0"
            cpu_hash = cpu_usage_matrix[timeval]
            if progname in cpu_hash:
                cpu_val = str(cpu_hash[progname])
            cpuvals.append(cpu_val)

            
        mem_usage_ofh.write("\t".join(memvals) + "\n")
        cpu_usage_ofh.write("\t".join(cpuvals) + "\n")

    mem_usage_ofh.close()
    cpu_usage_ofh.close()

    print("Done. See output files: " + str([mem_usage_output_filename, cpu_usage_output_filename]))
    
        
    
        

    

def parse_collectl_dat(collectl_dat_file, cpu_usage_matrix, memory_usage_matrix, times, prognames):

    prev_time = 0
    prognames_dict = dict()

    with open(collectl_dat_file) as f:
        for line in f:
            vals = re.split("\s+", line)
            day = vals[0]
            timeval_str = vals[1]
            timestruct = datetime.datetime.strptime("{} {}".format(day, timeval_str), "%Y%m%d %H:%M:%S")
            timeval_numeric = time.mktime(timestruct.timetuple())

            timeval_key_str = str(timeval_numeric)
            if timeval_numeric > prev_time:
                times.append(timeval_numeric)
                prev_time = timeval_numeric
            
            progname = os.path.basename(vals[19])
            cpu_usage_timehash = cpu_usage_matrix[timeval_key_str]
            if progname not in cpu_usage_timehash:
                cpu_usage_timehash[progname] = 1
            else:
                cpu_usage_timehash[progname] += 1
            
            
            memory = compute_GB(vals[8])
            progname_memory_hash = memory_usage_matrix[timeval_key_str]
            if progname not in progname_memory_hash:
                progname_memory_hash[progname] = memory
            else:
                progname_memory_hash[progname] += memory
            
            print("progname: {}, memoryG: {}".format(progname, memory))

            prognames_dict[progname] = 1

    progname_list = prognames_dict.keys()

    prognames += progname_list
        
    
    
def compute_GB(memory_val):

    r = re.search("^(\d+)([MKG])$", memory_val)
    numG = int(r.group(1))
    metric = r.group(2)

    if metric == 'M':
        numG /= 1e3
    elif metric == 'K':
        numG /= 1e6

    return numG 
    


 
####################
 
if __name__ == "__main__":
    main()
