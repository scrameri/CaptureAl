#!/bin/bash

## Usage: extract.fasta.seqs.sh -f <multi.fasta> -i <file with fasta sequence IDs (headers w/o ">")> 

## Needs:
# samtools

## Define arguments 
while getopts f:i: opts
do
        case "${opts}"
		in
 			f) fasta=${OPTARG};;
  			i) ids=${OPTARG};;
		esac
done

## Check arguments
if [ ! $fasta ] ; then echo "input FASTA not specified (-f option)" ; exit 0 ; fi
if [ ! $ids ] ; then echo "file with FASTA sequence IDs not specified (-i option)" ; exit 0 ; fi

if [ ! -f $fasta ] ; then echo "input FASTA <$fasta> not found" ; exit 0 ; fi
if [ ! -f $ids ] ; then echo "file with FASTA sequence IDs <$ids> not found" ; exit 0 ; fi

## Additional arguments
fbase=$(basename $fasta .fasta)

## Extract sequences
xargs samtools faidx $fasta < $ids > ${fbase}.extr.tmp

## Unwrap
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' ${fbase}.extr.tmp > ${fbase}.extr.fasta

## Clean up
rm -f ${fasta}.fai
rm -f ${fbase}.extr.tmp
