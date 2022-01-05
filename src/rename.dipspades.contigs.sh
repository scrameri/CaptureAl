#!/bin/bash

## Usage: rename.dipspades.contigs.sh -f <assembly folder> -s <sample file> -x <suffix> -t <number of threads>

## Needs: rename.dipspades.contigs.R

## Get arguments
while getopts f:s:x:t: opts
do
        case "${opts}"
        in
        	f) indir=${OPTARG};;
        	s) sfile=${OPTARG};;
		x) suffix=${OPTARG};;
		t) threads=${OPTARG};;
    	esac
done

## Check arguments
if [ ! $indir ] ; then echo "input folder not specified, stopping" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "input folder <${indir}> does not exist, stopping" ; exit 0 ; fi
if [ ! $sfile ] ; then echo "sample file not specified, stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <${sfile}> does not exist, stopping." ; exit 0 ; fi
if [ ! $threads ] ; then echo "number of threads not specified, setting to 8." ; threads=8 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping." ; exit 0 ; fi


## Define renaming function
doRename()
{
	id=$1
	echo $id
	rename.dipspades.contigs.R ${indir}/${id}${suffix} ${suffix} 1> /dev/null 2>> rename.err
}


## Execute renaming function
export -f doRename
export indir=$indir
export suffix=$suffix
if [ -f rename.err ] ; then /bin/rm -f rename.err ; fi
cat $sfile | parallel -j $threads doRename


## Finish
echo
echo "All samples processed!"
echo
