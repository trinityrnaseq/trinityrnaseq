#!/usr/bin/env python

import sys, os, string, subprocess

#aliasing the filenames using the labels


def run_command(command):
    print "Running command: " + command
    
    err_capture_file = open("my.stderr", 'w') # writing stderr to a file
    cmd_run = subprocess.Popen(args=command, shell=True, stderr=err_capture_file, stdout=sys.stdout)
    err = cmd_run.wait() # get exit code from command execution
    err_capture_file.close()

    if err:
        # report the error messages we captured, and exit non-zero
        sys.stderr.write("Error, cmd: " + command + " died with ret: " + `err`)
        for line in open(err_capture_file):
            sys.stderr.write(line)
        sys.exit(err)
    return

label_list = []  # symlink files to the labels
for i in range(1, len(sys.argv), 2):
    filename=sys.argv[i]
    label= sys.argv[i+1]
    cmd= "ln -sf " + filename + " " + label
    label_list.append(label)
    run_command(cmd)


# run the abundance estimation script

cmd = os.path.dirname(sys.argv[0]) + "/ctat_trinity_tool_wrapper.py " + " util/abundance_estimates_to_matrix.pl --gene_trans_map none --est_method RSEM --cross_sample_norm TMM " + " ".join(label_list)

run_command(cmd)

sys.exit(0)

