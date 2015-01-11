#!/bin/bash

cmd="./run_butterfly_nonJar.sh -N 10000 -L 300 -F 300 -C sample_data/RawComps.0/comp0 --stderr -V 20 "

eval $cmd

if [ "$?" -ne "0" ]; then
    echo "Error, command failed: " $cmd
    exit 1
fi

