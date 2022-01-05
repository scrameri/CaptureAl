#!/bin/bash

## Usage: create.multifastas.parallel.sh -l <locus file> -d <consensus dir paths> -p <prefix in consensus fasta relative to locus file> -s <suffix in consensus fasta relative to locus file> -t <number of threads>

## Define arguments
while getopts l:d:p:s:t: opts
do
        case "${opts}"
		in
				l) lfile=${OPTARG};;
				d) paths=${OPTARG};;
				p) prefix=${OPTARG};;
				s) suffix=${OPTARG};;
				t) threads=${OPTARG};;
		esac
done

## Check parameters and set defaults
if [ ! $lfile ] ; then echo "locus file not provided (-l option)" ; exit 0 ; fi
if [ ! -f $lfile ] ; then echo "locus file <$lfile> does not exist" ; exit 0 ; fi

if [ ! $paths ] ; then echo "file with paths to consensus directories (-d option) not specified" ; exit 0 ; fi
if [ ! -f $paths ] ; then echo "file with paths to consensus directories <$paths> does not exist" ; exit 0 ; fi

if [ ! $prefix ] ; then echo "consensus header prefix (-p option) not specifed, setting to ''" ; prefix='' ; fi
if [ ! $suffix ] ; then echo "consensus header suffix (-s option) not specified, setting to ''" ; suffix='' ; fi

if [ ! $threads ] ; then threads=4 ; echo "using 4 threads" ; else echo "using $threads threads" ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30" ; exit 0 ; fi

## Export arguments
export lfile=$lfile
export paths=$paths
export prefix=$prefix
export suffix=$suffix
nloc=$(wc -l < $lfile)
npaths=$(wc -l < $paths)

## Additional arguments
#export emptystring="--------------------------------------------------" # encodes a missing contig / sequence in output multifasta
export format=".fasta"                                                      # extension in input dir
export msuffix=".cons${format}"                                              # suffix of output multifasta

## Create output directory
export outdir="multifasta.cons.${nloc}.${npaths}" 
if [ -d $outdir ] ; then mv ${outdir} ${outdir}.bak ; echo "${outdir} exists, renaming it to ${outdir}.bak" ; fi
mkdir ${outdir}


## define multifasta function
doMultifasta()
{
	locus=$1
	echo ${locus}
	if [ -f ${locus}{msuffix} ] ; then /bin/rm -f ${locus}${msuffix} ; fi
	
	for dir in $(cat ../${paths})
	do
		fasta=$(ls -1 ../${dir}/${prefix}${locus}${suffix}${format})

		if [ ! -f $fasta ]
		then 
			echo "$fasta does not exist, stopping" >> ../multiconsensus.err
			exit 0
		else
			echo ">${dir}" >> ${locus}${msuffix}
			grep -v '^>' ${fasta} >> ${locus}${msuffix}
		fi
	done
}


## Execute function
export -f doMultifasta
cd ${outdir}
cat ../${lfile} | parallel -j $threads doMultifasta
cd ../


## Finish
echo ""
echo "All loci processed."
echo ""
