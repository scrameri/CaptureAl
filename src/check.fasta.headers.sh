#!/bin/bash

## Usage: check.fasta.headers <folder>

## Get argument
folder=$1

## Additional arguments
suffix=".fasta"

## Check argument
if [ ! $folder ] ; then echo "input folder not specified (1st and only argument)" ; exit 0 ; fi
if [ ! -d $folder ] ; then echo "input folder <$folder> does not exist" ; exit 0 ; fi

## Count number uf unique sequence names
for f in ${folder}/*${suffix} ; do grep '^>' ${f} ; done | sort -u > ${folder}${suffix}.headers

## Check identity of sequence names in each alignment (with reference to first alignment)
if [ -f ${folder}${suffix}.check ] ; then /bin/rm -f ${folder}${suffix}.check ; fi

test1=$(grep '^>' $(ls -1d ${folder}/*${suffix} | head -n1))

for f in ${folder}/*${suffix}
do
	test2=$(grep '^>' $f)
	if [ "$test1" == "$test2" ]
	then 
		echo -e "${f}\tpassed" >> ${folder}${suffix}.check
	else
		echo -e "${f}\tfailed" >> ${folder}${suffix}.check
	fi
done

## Report results
nfiles=$(ls -1d ${folder}/*${suffix} | wc -l)
nheaders=$(wc -l < ${folder}${suffix}.headers)

echo
echo "directory name:                           $folder"
echo "number of fasta files with suffix ${suffix}: ${nfiles} (number of loci)"
echo
echo "number of unique fasta headers:           ${nheaders} (number of individuals)"
echo
echo "number of passed / failed fasta files:"
awk '{print $2}' ${folder}${suffix}.check | sort | uniq -c
echo
