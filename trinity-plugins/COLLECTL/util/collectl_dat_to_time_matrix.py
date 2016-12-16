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

    parser.add_argument("--debug", required=False, action="store_true", default=False, help="debug mode")

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)      

    
    cpu_usage_matrix = collections.defaultdict(dict)
    memory_usage_matrix = collections.defaultdict(dict)
    times = list()
    prognames = list()

    parse_collectl_dat(args.dat_file, cpu_usage_matrix, memory_usage_matrix, times, prognames)

    
    print("\t" + "\t".join(prognames))
    for timeval in times:
        vals = list()
        timeval = str(timeval)
        vals.append(timeval)
        for progname in prognames:
            memory_hash = memory_usage_matrix[timeval]
            memory_val = "0"
            if progname in memory_hash:
                memory_val = str(memory_hash[progname])
            vals.append(memory_val)
        print("\t".join(vals))
    
             

    

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
