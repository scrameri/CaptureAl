#!/bin/bash

## Usage: get.final.loci.sh <folder with filtered alignments>

## Get arguments
trimmeddir=$1

## Check arguments
if [ ! $trimmeddir ] ; then echo "folder with final alignments not specified" ; exit 0 ; fi
if [ ! -d $trimmeddir ] ; then echo "folder with final alignments <$trimmeddir> not found" ; exit 0 ; fi

## Get loci
cd $trimmeddir
ls -1 *.fasta > ../final.loci.names
cd ../

## Finish
nloc=$(wc -l < final.loci.names)

echo
echo "final number of loci: $nloc"
echo
