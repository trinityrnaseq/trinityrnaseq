#!/usr/bin/env bash

docker build -t minimap2 .
docker tag minimap2 "$DOCKER"/minimap2:2.17
docker push "$DOCKER"/minimap2:2.17
