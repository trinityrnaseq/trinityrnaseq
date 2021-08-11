#!/bin/bash

set -ex

datadir="trinity_ext_sample_data_singularity"

if [ ! -e "$datadir" ]; then
    git clone https://github.com/trinityrnaseq/trinity_ext_sample_data.git ${datadir}
fi


version=`cat VERSION.txt`
simg="trinityrnaseq.v${version}.simg"

singularity exec -e -H /tmp -B `pwd`/${datadir}:/${datadir} \
       ${simg} \
       bash -c "cd /${datadir} && make"

