#!/bin/bash

## Usage: get.fasta.lengths.parallel.sh -d <folder with .fasta files> -t <number of threads>

## Needs: get.fasta.lengths.sh

## Get arguments
while getopts 'd:t:' flag; 
do 
  case "${flag}" in
    d) folder="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
  esac
done

## Check arguments
if [ ! $folder ] ; then echo "input folder not specified, stopping" ; exit 0 ; fi
if [ ! $threads ] ; then echo "number of threads not specified, setting to 4" ; threads=4 ; fi

if [ ! -d $folder ] ; then echo "input folder <$folder> does not exist, stopping" ; exit 0 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
export suffix=".fasta"
export fbase=$(basename $folder)

## Create output directory
if [ -d ${fbase}.lengths ] ; then echo "output directory <${fbase}.lengths> already exists, moving to ${fbase}.lengths.bak" ; mv ${fbase}.lengths ${fbase}.lengths.bak ; fi
mkdir ${fbase}.lengths

## Define assessment function
GetLength()
{
	fas=$1
	fasbase=$(basename $fas $suffix)
	echo ${fasbase}
	
	#get.fasta.length.R $fas ${fbase}.lengths/${fasbase}.lengths # slower
	get.fasta.lengths.sh -f $fas -o 'tab' -d ${fbase}.lengths/

}

## Execute assessment function
export -f GetLength
ls -1d ${folder}/*${suffix} | parallel -j $threads GetLength


## Make a table of all individuals and contig lengths
cat $(ls -1 ${fbase}.lengths/*.lengths | head -n1) | cut -f1 > ${fbase}.locus.lengths.tmp1
awk -F '\t' 'FNR == 1 && FNR != NR { sep = "\t" } { line[FNR] = line[FNR] sep $2 } END { for(i = 1; i <= FNR; ++i) { print line[i] } }' $(ls -1 ${fbase}.lengths/*.lengths) > ${fbase}.locus.lengths.tmp2
paste ${fbase}.locus.lengths.tmp1 ${fbase}.locus.lengths.tmp2 > ${fbase}.locus.lengths


## Clean up
/bin/rm -f ${fbase}.locus.lengths.tmp*


## Finish
echo
echo "All ${suffix} files processed!"
echo
