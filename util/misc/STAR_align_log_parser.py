#!/usr/bin/env python

import sys, os, re



def main():

    usage = "\n\n\tusage: {} STAR.Log.final.out [STAR.Log.final.out ...]\n\n".format(sys.argv[0])
    if len(sys.argv) < 2:
        exit(usage)


    
    tokens = list()
    
    processed_first = False
    files = sys.argv[1:]

    for log_filename in files:

        token_to_val = dict()

        with open(log_filename) as fh:
            for line in fh:
                line = line.rstrip()
                vals = line.split("|")
                if len(vals) != 2:
                    continue
                key = vals[0].strip()
                val = vals[1].strip()
                if not processed_first:
                    tokens.append(key)

                token_to_val[key] = val

        if not processed_first:
            # print header
            print("\t".join(tokens))
            processed_first = True

        
        vals = [token_to_val[x] for x in tokens]
        print("\t".join(vals))
    

    sys.exit(0)



if __name__=='__main__':
    main()
