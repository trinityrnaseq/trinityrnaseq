#!/bin/bash

if [ ! -e trinity_ext_sample_data ]; then
    git clone https://github.com/trinityrnaseq/trinity_ext_sample_data.git
fi


docker run --rm -v `pwd`/trinity_ext_sample_data:/trinity_ext_sample_data \
       trinityrnaseq/trinityrnaseq \
       bash -c 'cd /trinity_ext_sample_data && make'

