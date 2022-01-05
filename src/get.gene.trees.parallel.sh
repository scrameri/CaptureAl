#!/bin/bash

## Usage: get.gene.trees.parallel.sh -d <directory with alignments> -n <number of bootstrap replicates> -t <threads>

## Needs:
# raxmlHPC-SSE3
# remove.empty.alignment.taxa.R and R/ape library

## Define arguments
while getopts d:n:t: opts
do
	case "${opts}"
		in
			d) folder=${OPTARG};;
			n) boot=${OPTARG};;
			t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! $folder ] ; then echo "input directory with alignments not specified (-d option)" ; exit 0 ; fi
if [ ! $boot ] ; then echo "number of bootstrap replicates not specified (-n option), setting to 100" ; boot=100 ; fi

if [ ! -d $folder ] ; then echo "input directory with alignments <$folder> not found" ; exit 0 ; fi

if [ ! $threads ] ; then echo "number of threads not specified (-t option), setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping." ; exit 0 ; fi

## Additional arguments
export alnsuffix=".fasta"           # suffix of FASTA alignments in $folder
export outdir="${folder}.GENETREES" # output folder (subfolders for each locus will be created here)
export dir=$(pwd)                   # working directory
export ratemodel="GTRCAT"           # model of rate heterogeneity

## Export arguments
export folder=$folder
export boot=$boot

## Create output directory
if [ -d ${outdir} ] ; then echo "output directory <$outdir> already exists, moving to ${outdir}.bak" ; mv ${outdir} ${outdir}.bak ; fi
mkdir ${outdir}

## Define function
runRaxml()
{
	fasta=$1
	fname=$(basename $fasta ${alnsuffix})
	echo $fname
	
	# create locus output directory
	mkdir ${outdir}/${fname}
	cd ${outdir}/${fname}
		
	# delete missing taxa (these will cause RAxML abortion - for ASTRAL, the input gene trees can have missing taxa)
	remove.empty.alignment.taxa.R ${dir}/${fasta} ${fname}.full.fasta
	
	# delete missing sites (if alignments were subsetted, there might be empty alignment columns)
	~/bin/RAxML-8.2/raxmlHPC-SSE3 -f c -s ${fname}.full.fasta -m ${ratemodel} -n reduced > /dev/null 2> RAxMLread.err
	if [ -e ${fname}.unwrpd.fasta.reduced ] ; then input=${fname}.full.fasta.reduced ; else input=${fname}.full.fasta ; fi

	# run RAxML using 1 thread (sequential version) <ratemodel> rate heterogeneity model and <boot> bootstrap replicates
	randx=$RANDOM
	randp=$RANDOM
	~/bin/RAxML-8.2/raxmlHPC-SSE3 -f a -m ${ratemodel} -x ${randx} -p ${randp} -s ${input} -n cat.rax -N ${boot} > RAxML.err 2 > RAxML.log

	cd ../../
}

## Run RAxML on individual gene trees
export -f runRaxml
ls -1 ${folder}/*${alnsuffix} | parallel -j $threads runRaxml

## Collect gene trees
cat ${outdir}/*/RAxML_bipartitions.cat.rax > ${folder}.genetrees

## Finish
echo 
echo "All gene trees finished!"
echo
