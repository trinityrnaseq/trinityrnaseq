#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build trinityrnaseq.v${VERSION}.simg docker://trinityrnaseq/trinityrnaseq:$VERSION

singularity exec -e trinityrnaseq.v${VERSION}.simg Trinity --version


