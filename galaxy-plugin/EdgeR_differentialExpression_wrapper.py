import sys, os, subprocess

TRINITY_BASE_DIR = ""
if os.environ.has_key('TRINITY_HOME'):
    TRINITY_BASE_DIR = os.environ['TRINITY_HOME'];
else:
    sys.stderr.write("You must set the environmental variable TRINITY_BASE_DIR to the base installation directory of Trinity before running this");
    sys.exit(1)

usage= "usage: " + " $counts_matrix" + " $dispersion"

if len(sys.argv)<2:
    print "Require atleast two parameters"
else:
    print "All good- command going ahead"
print " "

def run_command(cmd):
    print "The command used: " + cmd
    pipe=subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    pipe.wait()
    ret= pipe.returncode
    if ret:
        print "command died: " + str(ret)
        print pipe.stderr.readlines()
        sys.exit(1)
    else:
        return
print " "

countmatrix= "counts_matrix"

cmd= "cp " + sys.argv[1] + " " + countmatrix 
run_command(cmd)

cmd= TRINITY_BASE_DIR + "/Analysis/DifferentialExpression/run_DE_analysis.pl "+ " --matrix "+ countmatrix + " --method edgeR " + " --output edgeR_results "+ " --dispersion " + sys.argv[2] + " --tar_gz_outdir"

run_command(cmd)

sys.exit(0)
