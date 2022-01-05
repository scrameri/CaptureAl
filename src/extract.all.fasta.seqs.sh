#!/bin/bash

## Usage: extract.all.fasta.seqs.sh <file.fasta>

## Needs: exonerate

## Get arguments
fasta=$1

## Check arguments
if [ ! $fasta ] ; then echo "input file not specified" ; exit 0 ; fi
if [ ! -f $fasta ] ; then echo "input file <$fasta> not found" ; exit 0 ; fi

## Additional arguments
refbase=$(basename $fasta .fasta)
dir=$(pwd)

# Index and extract reference sequences
if [ -d $refbase.seqs ] ; then echo "output directory <$refbase.seqs> exists, moving to ${refbase}.seqs.bak" ; mv ${refbase}.seqs ${refbase}.seqs.bak ; fi
mkdir ${refbase}.seqs
cd ${refbase}.seqs

echo
echo "extracting reference sequences..."

fastaindex ${dir}/${fasta} ${fbase}.idx

for locus in $(grep '^>' ${dir}/${fasta} | cut -f2 -d'>')
do
	echo $locus

	# fetch
	fastafetch -f ${dir}/${fasta} -i ${fbase}.idx -q ${locus} > ${locus}.fasta.tmp

	# unwrap
	cat ${locus}.fasta.tmp | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' > ${locus}.fasta
	/bin/rm -f ${locus}.fasta.tmp
done

echo "done"
echo

cd ../
