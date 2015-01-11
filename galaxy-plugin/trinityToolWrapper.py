#!/usr/bin/env python


# borrowed from: http://wiki.g2.bx.psu.edu/Future/Job%20Failure%20When%20stderr and modified for use with Trinity tools.

"""
Wrapper that execute a program and its arguments but reports standard error
messages only if the program exit status was not 0
Example: ./stderr_wrapper.py myprog arg1 -f arg2
"""

import sys, subprocess, os

assert sys.version_info[:2] >= ( 2, 4 )

TRINITY_BASE_DIR = ""
if os.environ.has_key('TRINITY_HOME'):
    TRINITY_BASE_DIR = os.environ['TRINITY_HOME'];
else:
    sys.write("You must set the environmental variable TRINITY_BASE_DIR to the base installation directory of Trinity before running this");
    sys.stderr.write("You must set the environmental variable TRINITY_BASE_DIR to the base installation directory of Trinity before running this");
    sys.exit(1)



# get bindir
bindir = sys.argv[0]
bindir = bindir.split("/")
if len(bindir) > 1:
    bindir.pop()
    bindir = "/".join(bindir)
else:
    bindir = "."


## add locations of tools to path setting.
TOOL_PATHS_FILE = bindir + "/__add_to_PATH_setting.txt";
for line in open(TOOL_PATHS_FILE):
    line = line.rstrip()
    os.environ['PATH'] += ":" + line


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

    args[0] = "".join([TRINITY_BASE_DIR, '/', args[0]]);

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
