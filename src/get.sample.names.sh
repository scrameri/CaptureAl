#!/bin/bash

## Usage: get.sample.names.sh <assembly folder> <sample suffix (will be removed in output)>

## Value: extracts all sample names in <assembly folder> and removes <sample suffix>

## Get arguments
assemdir=$1
samplesuffix=$2
oname="samples.txt"

## Extract
for sample in $(ls -1d ${assemdir}/*${samplesuffix})
do 
	basename $sample $samplesuffix
done > $oname
