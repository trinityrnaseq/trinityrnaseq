import sys, os, subprocess

TRINITY_BASE_DIR = ""
if os.environ.has_key('TRINITY_HOME'):
    TRINITY_BASE_DIR = os.environ['TRINITY_HOME'];
elif hasattr(os, 'symlink'): # symlink was implemented to always return false when it was not implemented in earlier versions of python.
    # 2017-09-26
    # Cicada Dennis added looking for the location of the Trinity program using the Unix "which" utility.
    # I tried using "command -v Trinity" but for some reason, I was getting a OS permission error with that.
    # I also found distutils.spawn.find_executable() which might work, but already implemented the below.
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
    sys.stderr.write("Either set TRINITY_HOME to the trinity base directory, " + \
        "or ensure that directory is in the PATH before running.")
    sys.exit(1)

usage= "usage: " + " $counts_matrix" + " $dispersion"

if len(sys.argv)<2:
    print "Require at least two parameters, " + usage
else:
    print "All good- command going ahead"
print " "

def run_command(cmd):
    # 2017-10-02
    # Cicada Dennis put the subprocess command in a try/except statement.
    # Errors were going undetected the way it was written previously.
    print "Running command:\n\t" + cmd
    try:
        pipe = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        cmd_stdout, cmd_stderr = pipe.communicate()
    except:
        msg = "ERROR while running command:\n\t" + cmd
        sys.stderr.write(msg)
        raise
    
    sys.stdout.write(cmd_stdout)
    ret = pipe.returncode
    if ret:
        print "command died: " + str(ret)
        sys.stderr.write(cmd_stderr)
        sys.exit(ret)
    else:
        # Any error output is written to stdout instead of stderr, since no error has occurred.
        sys.stderr.write(cmd_stderr) 
        return

print ""

countmatrix= "counts_matrix"

cmd= "cp " + sys.argv[1] + " " + countmatrix 
run_command(cmd)

cmd= TRINITY_BASE_DIR + "/Analysis/DifferentialExpression/run_DE_analysis.pl "+ " --matrix "+ countmatrix + " --method edgeR " + " --output edgeR_results "+ " --dispersion " + sys.argv[2] + " --tar_gz_outdir"

run_command(cmd)

sys.exit(0)
