#!/bin/bash

dirname=`dirname $0`

java -cp $dirname/bin:$dirname/lib/Jaligner.jar NWalign $*

exit $?

