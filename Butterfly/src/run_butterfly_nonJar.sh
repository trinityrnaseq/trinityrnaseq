#!/bin/bash

bindir=`dirname $0`

cmd="java -cp $bindir/bin:$bindir/lib/collections-generic-4.01.jar:$bindir/lib/java-getopt-1.0.13.jar:$bindir/lib/jung-algorithms-2.0.1.jar:$bindir/lib/jung-api-2.0.1.jar:$bindir/lib/jung-graph-impl-2.0.1.jar:$bindir/lib/Jaligner.jar TransAssembly_allProbPaths $*"

eval $cmd

exit $?


