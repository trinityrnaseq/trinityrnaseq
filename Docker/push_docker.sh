#!/bin/bash

VERSION=`cat VERSION.txt`

docker push trinityrnaseq/trinityrnaseq:${VERSION} 
docker push trinityrnaseq/trinityrnaseq:latest 

