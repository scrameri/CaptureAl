#!/bin/bash

## Score topology using gene trees and ASTRAL III
genetrees=$1
topology=$2
outname=$topology.score
mapfile=$3

if [ ! $mapfile ]
then
	echo "no mapfile provided"
	java -Xmx1G -jar ~/bin/ASTRAL_v5.6.3/astral.5.6.3.jar -i $genetrees -q $topology -o ${outname}.spectree -t 2 2> ${outname}.single.log
else
	echo "mapfile <$mapfile> provided"
	if [ ! -f $mapfile ] ; then echo "mapfile <$mapfile> not found, stopping." ; exit 0 ; fi
	java -Xmx1G -jar ~/bin/ASTRAL_v5.6.3/astral.5.6.3.jar -i $genetrees -q $topology -o ${outname}.spectree -a $mapfile -t 2 2> ${outname}.log
fi
