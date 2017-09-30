#!/usr/bin/env python

import os, re, sys, subprocess

trinity_package_dir = os.path.abspath(sys.argv[0] + "/../../../")
print("Trinity package dir: {}".format(trinity_package_dir))

trinity_package_name = os.path.basename(trinity_package_dir)
print("Trinity package name: {}".format(trinity_package_name))

destination_dir = "/usr/local/bin"
destination_package_dir = destination_dir + "/" + trinity_package_name

subprocess.check_call("rsync -av --exclude='.*' {} {}".format(trinity_package_dir, destination_package_dir), shell=True)

print("Trinity package installed at: {}".format(destination_package_dir))
print("\n\n\tFor convenience, set env var TRINITY_HOME={}".format(destination_package_dir))
print("\n\tSimply add:\n\n\texport TRINITY_HOME={}".format(destination_package_dir))
print("\tto your ~/.bashrc file.\n\n")
print("\tand run trinity via:  $TRINITY_HOME/Trinity\n\n\n")

sys.exit(0)




