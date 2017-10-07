#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t trinityrnaseq/trinityrnaseq:${VERSION} .
docker build -t trinityrnaseq/trinityrnaseq:latest .

