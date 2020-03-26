#!/bin/bash

## Usage: get.readcounts.parallel.sh -s <samples.txt> -p <prefix> -x <suffix> -t <threads>

## Authors: Simon Crameri (ETHZ), Stefan Zolelr (GDC) 

## Define arguments
while getopts s:p:x:t: opts
do
        case "${opts}"
		in
				s) sfile=${OPTARG};;
				p) prefix=${OPTARG};;
				x) suffix=${OPTARG};;
				t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! ${sfile} ] ; then echo "sample file not provided (-s option), stopping" ; exit 0 ; fi
if [ ! -f ${sfile} ] ; then echo "sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! ${prefix} ] ; then echo "sample prefix not provided (-p option), setting to ''" ; prefix="" ; fi
if [ ! ${suffix} ] ; then echo "sample suffix not provided (-x option), setting to '_R1.fastq.gz'" ; suffix="_R1.fastq.gz" ; fi
if [ ! ${threads} ] ; then echo "number of threads not provided (-t option), setting to 8" ; threads=8 ; fi
if [ ${threads} -gt 30 ] ; then echo "number of threads (-t option) must not exceed 30, stopping" ; exit 0 ; fi

## Check output file
if [ -f read.counts.txt ] ; then echo "read.counts.txt already exists, moving to read.counts.txt.bak" ; mv read.counts.txt read.counts.txt.bak ; fi

## Define function
get_readcounts()
{
	s=$1
	echo ${s}

	echo ${s} | tr "\n" "\t" > read.counts.${s}.tmp
	zcat ${prefix}${s}${suffix} |wc -l | tr "\n" "\t" |perl -n -e 'print $_ / 4' >> read.counts.${s}.tmp
	echo "" >> read.counts.${s}.tmp
	
}

## Execute function
export prefix=${prefix}
export suffix=${suffix}
export -f get_readcounts

cat ${sfile} | parallel -j ${threads} get_readcounts

## Combine
cat read.counts.*.tmp > read.counts.txt

## Clean up
rm -f read.counts.*.tmp


## Finish
echo
echo "All samples processed!"
echo
