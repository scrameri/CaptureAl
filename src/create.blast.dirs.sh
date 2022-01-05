#!/bin/bash

## Usage: create.sample.dirs.sh -s <sample file> -d <directory with assemblies for each sample> -x <sample suffix> -c <conserved path to contigs> -t <threads>

## Author: simon.crameri@env.ethz.ch, Apr 2019

## Define arguments
while getopts s:d:x:c:t: opts
do
        case "${opts}"
        in
        	s) sfile=${OPTARG};;
        	d) assemdir=${OPTARG};;
        	x) suffix=${OPTARG};;
        	c) contigspath=${OPTARG};;
        	t) threads=${OPTARG};;
        esac
done

## Check parameters and set defaults
if [ ! $sfile ] ; then echo "sample file not provided (-s option), stopping" ; exit 0 ; fi
if [ ! $assemdir ] ; then echo "directory with assemblies for each sample not provided (-d option), stopping" ; exit 0 ; fi
if [ ! $suffix ] ; then echo "sample suffix (string after {sample} in ${assemdir}) not provided (-x option), setting to ''" ; suffix='' ; fi
if [ ! $contigspath ] ; then echo "path from ${assemdir}/{sample}${suffix}/ to contigs fasta not provided (-c option), stopping" ; exit 0 ; fi

if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! -d $assemdir ] ; then echo "assembly directory <$assemdir> not found, stopping" ; exit 0 ; fi

if [ ! $threads ] ; then echo "number of threads not specified, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
export outdir="blast.contigs"    # name of output directory

## Export assembly directory and contigs path
# example1 for <myassemblies/sample1.dipspades/sample1.spades/dipspades/consensus_contigs.fasta>
# example2 for <myassemblies/sample1.spades/consensus_contigs.fasta>
export assemdir=$assemdir        # e.g. 'myassemblies', the directory containing subdirectories with a single contigspath
export suffix=$suffix            # e.g. '.dipspades' in example1 or '.spades' in example2
export contigspath=$contigspath  # e.g. '.spades/dipspades/consensus_contigs.fasta' in example1 or 'consensus_contigs.fasta' in example2

## Create output dir in current directory
if [ -d ${outdir} ]
then 
	echo "${outdir} already exists, moving it to ${outdir}.bak" 
	mv ${outdir} ${outdir}.bak
fi
mkdir ${outdir}

## Define function
makeContigsFolders()
{
	sample=$1
	echo ${sample}
	
	# looks for $contigspath assuming one of two possible folder architectures
	contigs="${assemdir}/${sample}${suffix}/${sample}${contigspath}"
	contigs2="${assemdir}/${sample}${suffix}/${contigspath}"
	
	# create output subdirectory
	mkdir ${outdir}/${sample}
	
	# copy contigs file
	if [ -f ${contigs} ]
	then
		cp ${contigs} ${outdir}/${sample}
	else
		if [ -f ${contigs2} ]
		then
			cp ${contigs2} ${outdir}/${sample}
		else
			echo "contigs file not found, check paths in ${outdir}/create.dirs.err"
			#echo ${sample}  >> ${outdir}/create.dirs.err
			echo ${contigs}  >> ${outdir}/create.dirs.err
			echo ${contigs2} >> ${outdir}/create.dirs.err
			echo ""          >> ${outdir}/create.dirs.err
		fi
	fi
}

## Execute function
export -f makeContigsFolders
cat ${sfile} | parallel -j $threads makeContigsFolders

## Finish
echo
echo "Done!"
echo
