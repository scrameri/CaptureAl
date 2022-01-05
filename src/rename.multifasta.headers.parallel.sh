#!/bin/bash

## Usage: rename-multifasta-headers-parallel.sh <folder with .fasta files> <nthreads> <pattern1> <pattern2> <stringsplit pattern>

## Get arguments
while getopts 'd:a:b:s:t:' flag; 
do
  case "${flag}" in
    d) dir="${OPTARG}";;
    a) pattern1="${OPTARG}";;
    b) pattern2="${OPTARG}";;
    s) splitchar="${OPTARG}";;
    t) threads="${OPTARG}";;
  esac
done

## Check arguments
if [ ! $dir ] ; then echo "input directory not specified (-d option)" ; exit 0 ; fi
if [ ! $pattern1 ] && [ ! $pattern2 ] ; then echo "no strings removed (-a option, -b option)" ; fi
if [ ! $pattern1 ] ; then pattern1='FALSE' ; fi
if [ ! $pattern2 ] ; then pattern2='FALSE' ; fi
if [ ! $splitchar ] ; then echo "split string not specified (-s option), no string splitting done" ; splitchar='FALSE' ; fi
if [ "$pattern1" == "FALSE" ] && [ "$pattern2" == "FALSE" ] && [ "$splitchar" == "FALSE" ] ; then exit 0 ; fi

if [ ! -d $dir ] ; then echo "input directory <$dir> does not exist" ; exit 0 ; fi

if [ ! $threads ] ; then echo "number of threads not specified (-t option), setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30" ; exit 0 ; fi

## Additional arguments
suffix=".fasta" # will look for files with this extension in $dir

## Export 
export pattern1=$pattern1
export pattern2=$pattern2
export splitchar=$splitchar

## Define
doRename()
{

	file=$1
	rename.fasta.headers.R $file $pattern1 $pattern2 $splitchar

}

## Execute
export -f doRename
ls -1 ${dir}/*${suffix} | parallel -j $threads doRename

## Finish
echo
echo "All $suffix files processed!"
echo