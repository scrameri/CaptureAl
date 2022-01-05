#!/bin/bash

## Usage: symlink.bams.sh -s <samples.txt> -Q <mapping quality threshold> -b <path to BAM files> -o <output directory>

# -s sample file
# -Q mapping quality threshold
# -b path to BAM files
# -o output directory

## Define arguments
while getopts s:Q:b:o: opts
do
        case "${opts}"
        in
       		s) sfile=${OPTARG};;
       		Q) Q=${OPTARG};;
			b) bampath=${OPTARG};;
			o) odir=${OPTARG};;
    	 esac
done


## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) required, stopping." ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping." ; exit 0 ; fi

if [ ! $Q ] ; then echo "mapping quality parameter (-Q option) not set, assuming -Q 10." ; Q=10 ; fi

if [ ! $bampath ] ; then echo "path to BAM file not provided (-b option), setting to <SAMPLE/SAMPLE.bwa-mem.mapped.Q${Q}.sorted.bam>." ; bampath="SAMPLE/SAMPLE.bwa-mem.mapped.Q${Q}.sorted.bam" ; fi

if [ ! $odir ] ; then echo "path to output directory (-o option) required, stopping." ; exit 0 ; fi


## Additional arguments

## Symlink BAM and BAI files
if [ ! -d $odir ] ; then mkdir $odir ; fi

cd $odir

for sample in $(cat ../$sfile)
do
	# set path to BAM file
	bam=$(echo "${bampath}" | sed -e "s/SAMPLE/$sample/g")
	
	# symlink
	ln -s ../${bam} .
	ln -s ../${bam}.bai .
done

cd ../