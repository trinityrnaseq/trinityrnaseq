#!/bin/bash

echo "## Checking plugin installations:"
echo

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
