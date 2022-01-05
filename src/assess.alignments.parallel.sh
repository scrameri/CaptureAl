#!/bin/bash

## Usage: assess.alignments.parallel.sh -f <folder with .fasta alignments> -t <number of threads>

## Needs:
# assess.alignment.R (needs R ape library)

## Define arguments
while getopts f:t: opts
do
	case "${opts}"
		in
			f) folder=${OPTARG};;
			t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! ${folder} ] ; then echo "folder with alignments (-f option) not specified, stopping." ; exit 0 ; fi
if [ ! -d ${folder} ] ; then echo "folder with alignments <${folder}> does not exist, stopping." ; exit 0 ; fi

if [ ! ${threads} ] ; then echo "number of threads (-t option) not specified, setting to 4" ; threads=4 ; fi
if [ ${threads} -gt 30 ] ; then echo "number of threads must not exceed 30" ; exit 0 ; fi

## Additional arguments
export alnsuffix=".all.aln.etr.itr.fasta"
export outdir="${folder}.assess"
ofile="${folder}.assess.txt"

## Create output dirs and move any existing dirs or files
if [ -d ${outdir} ] ; then mv ${outdir} ${outdir}.bak ; fi
mkdir ${outdir}

## Define assessment function
doAssess()
{
	aln=$1
	alnbase=$(basename $aln ${alnsuffix})
	echo $alnbase
	
	assess.alignment.R $aln ${outdir}/${alnbase}.assess
}

## Execute assessment function
export -f doAssess
ls -1d ${folder}/*${alnsuffix} | parallel -j $threads doAssess

## Concatenate
echo "Concatenating results..."
head -n1 $(ls -1d ${outdir}/* | head -n1) > ${ofile}
for file in $(ls -1d ${outdir}/*)
do 
tail -n+2 $file >> ${ofile}
done

## Clean up
/bin/rm -rf ${outdir}

## Finish
echo
echo "All alignments assessed!"
echo
