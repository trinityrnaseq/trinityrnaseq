#!/bin/bash

if [ ! -e trinity_ext_sample_data ]; then
    git clone https://github.com/trinityrnaseq/trinity_ext_sample_data.git
fi


VERSION=`cat VERSION.txt`

docker run --rm -v `pwd`/trinity_ext_sample_data:/trinity_ext_sample_data \
       trinityrnaseq/trinityrnaseq:$VERSION \
       bash -c 'cd /trinity_ext_sample_data && make'

