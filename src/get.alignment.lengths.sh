#!/bin/bash

## Usage: get.alignment.lengths.sh <folder with alignments>

## Get arguments
algn=$1

## Additional arguments
alnsuffix=".fasta"

## Remove $algn.names.tmp | $algn.length.tmp | $algn.lengths if present
if [ -f $algn.length.tmp ] ; then /bin/rm -f $algn.length.tmp ; fi
if [ -f $algn.names.tmp ] ; then /bin/rm -f $algn.names.tmp ; fi
if [ -f $algn.lengths ] ; then mv ${algn}.lengths ${algn}.lengths.bak ; fi

## Loop through alignments
ls -1 $algn/*$alnsuffix > $algn.names.tmp

for fasta in $(cat $algn.names.tmp)
do 
	echo $(basename $fasta .fasta)
	awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $fasta | head -n 2 | tail -n 1 >> $algn.length.tmp
done

## Paste .fasta name and alignment lengths
paste $algn.names.tmp $algn.length.tmp > ${algn}.lengths

## Clean up
/bin/rm -f $algn.names.tmp
/bin/rm -f $algn.length.tmp
