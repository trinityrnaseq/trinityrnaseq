#!/bin/bash

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo 'Performing Unit Tests of Build'
echo ' '
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

if [ -e "trinity-plugins/jellyfish/jellyfish" ]
then
	echo "JellyFish:               has been Installed Properly"
else
	echo "JellyFish                Installation appears to have FAILED"
fi

if [ -e "Inchworm/bin/inchworm" ]
then
	echo "Inchworm:                has been Installed Properly"
else
	echo "Inchworm Installation appears to have FAILED"
fi
if [ -e "Chrysalis/Chrysalis" ]
then
	echo "Chrysalis:               has been Installed Properly"
else
	echo "Chrysalis Installation appears to have FAILED"
fi
if [ -e "Chrysalis/QuantifyGraph" ]
then
	echo "QuantifyGraph:           has been Installed Properly"
else
	echo "QuantifyGraph Installation appears to have FAILED"
fi

if [ -e "Chrysalis/GraphFromFasta" ]
then
	echo "GraphFromFasta:          has been Installed Properly"
else
	echo "GraphFromFasta Installation appears to have FAILED"
fi

if [ -e "Chrysalis/ReadsToTranscripts" ]
then
	echo "ReadsToTranscripts:      has been Installed Properly"
else
	echo "ReadsToTranscripts Installation appears to have FAILED"
fi

if [ -e "trinity-plugins/fastool/fastool" ]
then
	echo "fastool:                 has been Installed Properly"
else
	echo "fastool Installation appears to have FAILED"
fi
if [ -e "trinity-plugins/parafly/bin/ParaFly" ]
then
	echo "parafly:                 has been Installed Properly"
else
	echo "parafly Installation appears to have FAILED"
fi

if [ -e "trinity-plugins/BIN/samtools" ]
then
    echo "samtools-0.1.19          has been Installed Properly"
else
    echo "samtools-0.1.19 Installation appears to have FAILED"
fi
