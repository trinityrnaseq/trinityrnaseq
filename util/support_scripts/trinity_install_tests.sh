#!/bin/bash

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo 'Performing Unit Tests of Build'
echo ' '
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"


if [ -e "Inchworm/bin/inchworm" ]
then
	echo "Inchworm:                has been Installed Properly"
else
	echo "Inchworm Installation appears to have FAILED"
fi
if [ -e "Chrysalis/bin/Chrysalis" ]
then
	echo "Chrysalis:               has been Installed Properly"
else
	echo "Chrysalis Installation appears to have FAILED"
fi
if [ -e "Chrysalis/bin/QuantifyGraph" ]
then
	echo "QuantifyGraph:           has been Installed Properly"
else
	echo "QuantifyGraph Installation appears to have FAILED"
fi

if [ -e "Chrysalis/bin/GraphFromFasta" ]
then
	echo "GraphFromFasta:          has been Installed Properly"
else
	echo "GraphFromFasta Installation appears to have FAILED"
fi

if [ -e "Chrysalis/bin/ReadsToTranscripts" ]
then
	echo "ReadsToTranscripts:      has been Installed Properly"
else
	echo "ReadsToTranscripts Installation appears to have FAILED"
fi


if [ -e "trinity-plugins/BIN/ParaFly" ]
then
	echo "parafly:                 has been Installed Properly"
else
	echo "parafly Installation appears to have FAILED"
fi

