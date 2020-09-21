#!/usr/bin/env bash

docker build -t trinityrnaseq-wdl .
docker tag trinityrnaseq-wdl "$DOCKER"/trinityrnaseq-wdl:1.0.0
docker push "$DOCKER"/trinityrnaseq-wdl:1.0.0
