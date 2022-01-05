#!/bin/bash

## Usage: subset.alignments.sh -s <sample file> -d <dir with alignments> -o <output directory> -g <FLAG whether to remove alignment columns entirely consisting of gaps> -t <number of threads> 

## Needs:
# subset.alignment.R (needs R ape library)

## Define arguments
ignoregaps="false"
while getopts s:d:o:t:g opts
do
	case "${opts}"
		in
			s) sfile=${OPTARG};;
			d) indir=${OPTARG};;			
			o) outdir=${OPTARG};;
			g) ignoregaps="true";;
			t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! ${sfile} ] ; then echo "sample file (-s option) not specified" ; exit 0 ; fi
if [ ! -f ${sfile} ] ; then echo "sample file  <${sfile}> does not exist" ; exit 0 ; fi

if [ ! ${indir} ] ; then echo "input directory (-d option) not specified" ; exit 0 ; fi
if [ ! -d ${indir} ] ; then echo "input directory <${indir}> does not exist" ; exit 0 ; fi

if [ ! ${outdir} ] ; then echo "output directory (-o option) not specified, setting to ${indir}.sub" ; outdir="${indir}.sub" ; fi

if [ ! ${threads} ] ; then echo "number of threads (-t option) not specified, setting to 4" ; threads=4 ; fi
if [ ${threads} -gt 30 ] ; then echo "number of threads must not exceed 30" ; exit 0 ; fi

## Additional arguments
export alnsuffix=".fasta"

## Create output dir
if [ -d ${outdir} ] ; then echo "${outdir} already exists, moving ot ${outdir}.bak" ; mv ${outdir} ${outdir}.bak ; fi
mkdir ${outdir}

## Export arguments
export sfile=$sfile
export indir=$indir
export outdir=$outdir
export ignoregaps=$ignoregaps

## Define function
do_subset()
{
	aln=$1
	alnbase=$(basename $aln ${alnsuffix})
	echo $alnbase
	
	#oname="${alnbase}.sub${alnsuffix}"
	oname="${alnbase}${alnsuffix}"

	args="$aln $sfile $ignoregaps $oname $outdir"
	subset.alignment.R $args
}

## Execute function
export -f do_subset
ls -1d ${indir}/*${alnsuffix} | parallel -j ${threads} do_subset

## Finish
echo 
echo "All alignments processed!"
echo
