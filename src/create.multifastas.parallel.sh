#!/bin/bash

## Usage: create.multifastas.parallel.sh -s <sample file> -l <locus file> -d <directory with best-scoring contigs> -t <number of threads>

## Define arguments
while getopts s:l:d:t: opts
do
        case "${opts}"
		in
				s) sfile=${OPTARG};;
				l) lfile=${OPTARG};;
				d) indir=${OPTARG};;	
				t) threads=${OPTARG};;
		esac
done

## Check parameters and set defaults
if [ ! $indir ] ; then echo "directory with best-scoring contigs not provided (-d option), stopping" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "directory <$indir> does not exist, stopping" ; exit 0 ; fi
if [ ! $lfile ] ; then echo "locus file not provided (-l option), stopping" ; exit 0 ; fi
if [ ! -f $lfile ] ; then echo "locus file <$lfile> does not exist, stopping" ; exit 0 ; fi
if [ ! $sfile ] ; then echo "samples file not provided (-s option), stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "samples file <$sfile> does not exist, stopping" ; exit 0 ; fi
if [ ! $threads ] ; then echo threads=15 ; "using 15 threads" ; else echo "using $threads threads" ; fi

## Export arguments
export sfile=$sfile
export lfile=$lfile
export indir=$indir

nind=$(wc -l < $sfile)
nloc=$(wc -l < $lfile)


## Additional arguments
#export emptystring="--------------------------------------------------" # encodes a missing contig / sequence in output multifasta
export msuffix=".all.fasta"                                              # suffix of output multifasta
export fsuffix=".bestScore.fasta"                                        # suffix of input fasta
export dir=$(pwd)


## Create output directory
export outdir="multifasta.${nind}.${nloc}" 
if [ -e $outdir ] ; then mv ${outdir} ${outdir}.bak ; echo "${outdir} exists, renaming it to ${outdir}.bak" ; fi
mkdir ${outdir}


## define multifasta function
doMultifasta()
{
	locus=$1
	echo ${locus}
	if [ -f ${locus}{msuffix} ] ; then /bin/rm -f ${locus}${msuffix} ; fi
	
	for sample in $(cat ${dir}/${sfile})
	do

		fasta=$(ls -1 ${dir}/${indir}/${sample}/${sample}.${locus}${fsuffix})

		if [ ! -f $fasta ]
		then 
			echo "$fasta does not exist, stopping" >> ${dir}/multifastas.err
			exit 0
		else
			cat ${fasta} >> ${locus}${msuffix}
		fi
	done
}


## Execute function
export -f doMultifasta
cd ${outdir}
cat ${dir}/${lfile} | parallel -j $threads doMultifasta
cd ${dir}


## Finish
echo ""
echo "All loci processed."
echo ""
