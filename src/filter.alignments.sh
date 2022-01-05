#!/bin/bash

## Usage: filter.alignments.sh -d <dir with alignments> -f <filterfile with alignment file names>

## Get arguments
while getopts 'd:f:' flag; 
do
  case "${flag}" in
    d) dir="${OPTARG}" ;;
    f) filterfile="${OPTARG}" ;;
  esac
done

## Check arguments
if [ ! $dir ] ; then echo "input folder with alignments not specified (-d option)" ; exit 0 ; fi
if [ ! $filterfile ] ; then echo "file with alignments to be filtered not specified (-f option)" ; exit 0 ; fi

if [ ! -d $dir ] ; then echo "input folder <$dir> not found" ; exit 0 ; fi
if [ ! -f $filterfile ] ; then echo "filter file <$filterfile> not found" ; exit 0 ; fi

## Additional arguments
nremoved=$(wc -l < $filterfile )
odirbase="loci_rm_"

## Reverse any previous filterings
nrmdir=$(ls -1d ${dir}/${odirbase}* 2> /dev/null | wc -l)
if [ $nrmdir -gt 0 ]
then
	for rmdir in $(ls -1d ${dir}/${odirbase}*)
	do
		mv ${rmdir}/* ${dir}/ 2> /dev/null
		/bin/rm -r ${rmdir}
	done
fi

## Create directory with alignments to be removed
if [ $nremoved -gt 0 ]
then
	filterbase=$(basename ${filterfile} .txt)
	mkdir ${dir}/${odirbase}${filterbase}
	for aln in $(cat ${filterfile})
	do 
		alnbase=$(echo ${aln} | awk '{ gsub(/.fasta$/, ""); print }')
		mv ${dir}/${alnbase}* ${dir}/${odirbase}${filterbase}/
	done
	echo
	echo "filtered $nremoved alignments!"
	echo
else
	echo
	echo "no alignments to be filtered!"
	echo
fi