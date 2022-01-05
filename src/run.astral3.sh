#!/bin/bash

## Usage run.astral3.sh -g <.genetrees> -s <species .mapfile>

## Define arguments
while getopts g:s: opts
do
	case "${opts}"
		in
			g) genetrees=${OPTARG};;
			s) mapspec=${OPTARG};;
		esac
done

## Check arguments
if [ ! $genetrees ] ; then echo "genetree file not specified (-g option)" ; exit 0 ; fi
if [ ! -f $genetrees ] ; then echo "genetree file <$genetrees> not found" ; exit 0 ; fi

if [ ! $mapspec ]
then 
	echo "no species mapfile provided, treating all individuals as different species" 
	single=1
else
	if [ ! -f $mapspec ] ; then echo "species mapfile <$mapspec> not found" ; exit 0 ; fi
	echo "running in multi-individual mode"
	single=0
fi

## Additional arguments
outname=$(basename $genetrees .genetrees)

## Run ASTRAL on multi-individual dataset using full annotation (-t 2)
if [ $single -eq 0 ]
then 	
	java -Xmx4G -jar ~/bin/ASTRAL_v5.6.3/astral.5.6.3.jar -i $genetrees -a ${mapspec} -o ${outname}.spectree -t 2 2> ${outname}.log
fi 

## Run ASTRAL on single-individual dataset using full annotation (-t 2)
if [ $single -eq 1 ]
then
	java -Xmx4G -jar ~/bin/ASTRAL_v5.6.3/astral.5.6.3.jar -i $genetrees -o ${outname}.single.spectree -t 2 2> ${outname}.single.log
fi

## Finish
echo
echo "Done ASTRAL analysis!"
echo
