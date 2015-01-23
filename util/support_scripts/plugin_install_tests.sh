#!/bin/bash

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo 'Performing Unit Tests of Build'
echo ' '
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

if [ -e "trinity-plugins/slclust/bin/slclust" ]
then
	echo "slclust:                 has been Installed Properly"
else
	echo "slclust Installation appears to have FAILED"
fi
if [ -e "trinity-plugins/collectl/bin/collectl" ]
then
	echo "collectl:                has been Installed Properly"
else
	echo "collectl Installation appears to have FAILED"
fi
if [ -e "trinity-plugins/rsem/rsem-calculate-expression" ]
then
    echo "rsem:                    has been Installed Properly"
else
    echo "rsem:    Installation appears to have FAILED"
fi

