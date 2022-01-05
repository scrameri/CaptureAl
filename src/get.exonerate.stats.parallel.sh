#!/bin/bash

## Usage: get.exonerate.stats.parallel.sh -s <samplefile> -d <path to dir with exonerate results> -t <number of threads>

## Define arguments
while getopts s:d:t: opts
do
        case "${opts}"
        in
                s) sfile=${OPTARG};;
                d) dir=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping!" ; exit 0 ; fi
if [ ! $dir ] ; then echo "directory with exonerate results (-d option) not specified, stopping!" ; exit 0 ; fi

if [ ! -f $sfile ] ; then echo "file <${sfile}> does not exist, stopping!" ; exit 0 ; fi
if [ ! -d $dir ] ; then echo "directory <${dir}> does not exist, stopping!" ; exit 0 ; fi

if [ ! $threads ] ; then echo "number of threads not specified, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping!" ; exit 0 ; fi

## Additional arguments
suffix="" # sample suffix in $dir

## Export
export dir=$dir

## Define function
doGetStats()
{
	sample=$1
	echo $sample
	
	get.exonerate.stats.R ${dir}/${sample}${suffix}
}

## Execute function
export -f doGetStats
cat $sfile | parallel -j $threads doGetStats

## Finish
echo
echo "All samples processed!"
echo

