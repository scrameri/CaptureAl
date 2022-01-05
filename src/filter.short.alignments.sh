#!/bin/bash

## Usage: filter.short.alignments.sh -d <folder with alignments> -f <file with alignment name (1st column) and alignment length (2nd column)> -m <minimum alignment length>

## Needs: 
# awk

## Get arguments
while getopts 'd:f:m:' flag; 
do 
  case "${flag}" 
  in
    d) dir="${OPTARG}" ;;
    f) lenfile="${OPTARG}" ;;
    m) minalnlen="${OPTARG}" ;;
  esac
done

## Check arguments
if [ ! $dir ] ; then echo "folder with alignments not specified (-d option)" ; exit 0 ; fi
if [ ! $lenfile ] ; then echo "file with alignment lenghts not specified (-f option)" ; exit 0 ; fi

if [ ! -d $dir ] ; then echo "folder with alignments <$dir> not found" ; exit 0 ; fi
if [ ! -f $lenfile ] ; then echo "file with alignment lenghts not found" ; exit 0 ; fi

if [ ! $minalnlen ] ; then echo "minimum alignment length not specified (-m option)" ; exit 0 ; fi
if [[ ! $minalnlen =~ ^[0-9]+$ ]] ; then echo "minimum alignment length <$minalnlen> is not a positive integer" ; exit 0 ; fi

## Additional arguments
odirbase="loci_rm_shorter_"

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

## Find short alignments
cat ${lenfile} | awk -v thresh=$minalnlen '($2 < thresh)' > short.aln.rm
nremoved=$(wc -l < short.aln.rm)

## Find shortest alignment
shortest=$(cat ${lenfile} | cut -f2 | sort -n | uniq | head -n1)
echo "shortest alignment is ${shortest}"

## Exit if shortest is >= minalnlen
if [ ! $nremoved -gt 0 ] ; then echo "no alignments are shorter than ${minalnlen}" ; /bin/rm -f short.aln.rm ; exit 0 ; fi

## Remove short alignments
mkdir ${dir}/${odirbase}${minalnlen}

for aln in $(cut -f1 short.aln.rm)
do
	mv $aln ${dir}/loci_rm_shorter_${minalnlen}
done

## Finish
echo
echo "removed $nremoved alignments shorter than ${minalnlen}"
echo 

## Clean up
/bin/rm -f short.aln.rm
