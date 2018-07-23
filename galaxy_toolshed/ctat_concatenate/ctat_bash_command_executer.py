#!/usr/bin/env python

import sys, subprocess, os

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

def __main__():
    # Get command-line arguments
    args = sys.argv
    # Remove name of calling program, i.e. ./stderr_wrapper.py
    args.pop(0)

    # If there are no arguments left, we're done
    if len(args) == 0:
        return

    # If one needs to silence stdout 
    #args.append( ">" )
    #args.append( "/dev/null" )

    cmdline = " ".join(args)
    
    
    try:
        # Run program
        err_capture = open("stderr.txt", 'w')
        proc = subprocess.Popen( args=cmdline, shell=True, stderr=err_capture, stdout=sys.stdout )
        returncode = proc.wait()
        err_capture.close()

        
        if returncode != 0:
            raise Exception

    except Exception:
        # Running Grinder failed: write error message to stderr
        err_text = open("stderr.txt").readlines()
        stop_err( "ERROR:\n" + "\n".join(err_text))


if __name__ == "__main__": __main__()
