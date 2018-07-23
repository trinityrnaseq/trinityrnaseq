#!/usr/bin/env python


# borrowed from: http://wiki.g2.bx.psu.edu/Future/Job%20Failure%20When%20stderr and modified for use with Trinity tools.

"""
Wrapper that execute a program and its arguments but reports standard error
messages only if the program exit status was not 0
Example: ./stderr_wrapper.py myprog arg1 -f arg2
"""

import sys, subprocess, os
print sys.version

assert sys.version_info[:2] >= ( 2, 4 )

TRINITY_BASE_DIR = ""
if os.environ.has_key('TRINITY_HOME'):
    TRINITY_BASE_DIR = os.environ['TRINITY_HOME'];
elif hasattr(os, 'symlink'): # symlink was implemented to always return false when it was not implemented in earlier versions.
    # 2017-09-26
    # Cicada Dennis added looking for the location of the Trinity program using the Unix "which" utility.
    # I tried using "command -v Trinity" but for some reason, I was getting a OS permission error with that.
    # I just found distutils.spawn.find_executable() which might work, but already implemented the below.
    try:
        pipe1 = subprocess.Popen(["which", "Trinity"], stdout=subprocess.PIPE)
    except:
        msg = "You must set the environmental variable TRINITY_HOME to the base installation directory of Trinity before running {:s}.\n".format(sys.argv[0])
        sys.stderr.write(msg)
        # t, v, tb = sys.exc_info()
        # raise t, v, tb
        # For some reason the above was giving a syntax error. 
        # A simple raise should reraise the existing exception.
        raise
    else:
        TrinityPath, err_info = pipe1.communicate()
        # FIX - probably should be checking err_info for errors...
        # Determine the TRINITY_BASE_DIR from output1.
        # If TrinityPath is a link, we need to dereference the link.
        TrinityPath = TrinityPath.rstrip() # Need to strip off a newline.
        # print "Trinity that was found is: {:s}".format(repr(TrinityPath))
        # print os.path.islink(TrinityPath)
        TrinityPath = os.path.abspath(TrinityPath)
        # msg = "The Absolute Trinity path that was found is: {:s}".format(TrinityPath)
        # print msg
        # print os.path.islink(TrinityPath)
        while os.path.islink(TrinityPath):
            # print "That path is a link."
            TrinityPath = os.path.join(os.path.dirname(TrinityPath),os.readlink(TrinityPath))
            # print "The new path is: {:s}".format(TrinityPath)
        # Take off the last part of the path (which is the Trinity command)
        TRINITY_BASE_DIR = "/".join(TrinityPath.split("/")[0:-1])
else:
    sys.stderr.write("Either set TRINITY_HOME to the trinity base directory, or ensure that directory is in the PATH before running.")
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
#TOOL_PATHS_FILE = bindir + "/__add_to_PATH_setting.txt";
#for line in open(TOOL_PATHS_FILE):
#    line = line.rstrip()
#    os.environ['PATH'] += ":" + line

# Add TrinityPath and its utils to the PATH environment variable.
# print "Initially the PATH env variable is:\n\t{:s}".format(os.environ['PATH'])
os.environ['PATH'] = os.environ['PATH'] + ":{:s}:{:s}".format(TRINITY_BASE_DIR,TRINITY_BASE_DIR+"/util")
# print "Now the PATH env variable is:\n\t{:s}".format(os.environ['PATH'])


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
    print "The TRINITY_BASE_DIR is:\n\t{:s}".format(TRINITY_BASE_DIR)
    print "The PATH env variable is:\n\t{:s}".format(os.environ['PATH'])

    args[0] = "".join([TRINITY_BASE_DIR, '/', args[0]]);

    cmdline = " ".join(args)

    print "The command being invoked is:\n\t{:s}".format(cmdline)
    

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
