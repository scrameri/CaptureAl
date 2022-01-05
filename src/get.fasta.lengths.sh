#!/bin/bash

## Usage: get.fasta.lengths.sh -f <fasta file> -o <optional: output format: 'fasta' or 'tab' (default)>

## Define input
while getopts f:o:d: opts
do
        case "${opts}"
        in
                f) fasta=${OPTARG};;
                o) format=${OPTARG};;
		d) outdir=${OPTARG};;
         esac
done

## Check input
if [ ! ${fasta} ] ; then echo "fasta file not specified (-f option), stopping" ; exit 0 ; fi
if [ ! -f ${fasta} ] ; then echo "fasta file <${fasta}> not found, stopping" ; exit 0 ; fi
if [ ! ${outdir} ] ; then outdir=$(pwd) ; fi
#if [ ! -d ${outdir} ] ; then echo "output directory <${outdir}> not found, stopping" ; exit 0 ; fi
if [ ! ${format} ] ; then format='tab' ; fi

arr='fasta tab'
if echo ${arr[@]} | grep -q -w ${format}
  then     
	 :
  else      
	echo "incorrect output format (-o option: specify as 'fasta' or 'tab')"
	exit 0
  fi

## Create output directory
fname=$(basename $fasta .fasta)
if [ ! -d ${outdir} ] ; then mkdir ${outdir} ; fi

## Fast way to get sequence lengths
# http://stackoverflow.com/questions/23992646/sequence-length-of-fasta-file
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $fasta > ${outdir}/${fname}.namelen.tmp

## Output
if [[ $format = 'fasta' ]]
then

	mv ${outdir}/${fname}.namelen.tmp ${outdir}/${fname}.lengths

else

	grep -v '>' ${outdir}/${fname}.namelen.tmp > ${outdir}/${fname}.len.tmp
	grep '>' ${outdir}/${fname}.namelen.tmp | cut -d ">" -f2 > ${outdir}/${fname}.name.tmp
	paste ${outdir}/${fname}.name.tmp ${outdir}/${fname}.len.tmp > ${outdir}/${fname}.lengths
	rm ${outdir}/${fname}.*.tmp

fi
