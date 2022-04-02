#!/bin/bash

## Usage: combine.contigs.parallel.sh -s <samplefile> -d <path to dir with assembly results> -e <path to dir with exonerate results> -a <minimum target alignment length> -c <minimum normalized alignment score> -t <number of threads>

## Define arguments
while getopts s:d:e:a:c:t: opts
do
        case "${opts}"
        in
                s) sfile=${OPTARG};;
                d) ass=${OPTARG};;
                e) dir=${OPTARG};;
                a) minaln=${OPTARG};;
                c) minnormscore=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping!" ; exit 0 ; fi
if [ ! $ass ] ; then echo "directory with assembly results (-d option) not specified, stopping!" ; exit 0 ; fi
if [ ! $dir ] ; then echo "directory with exonerate results (-e option) not specified, stopping!" ; exit 0 ; fi

if [ ! -f $sfile ] ; then echo "file <${sfile}> does not exist, stopping!" ; exit 0 ; fi
if [ ! -d $ass ] ; then echo "directory <${ass}> does not exist, stopping!" ; exit 0 ; fi
if [ ! -d $dir ] ; then echo "directory <${dir}> does not exist, stopping!" ; exit 0 ; fi

if [ ! $minaln ] ; then echo "minimum target alignment length set to minaln=80" ; minaln=80 ; fi
if [ ! $minnormscore ] ; then echo "minimum normalized alignment score set to minnormscore=2" ; minnormscore=2 ; fi

if [ ! $threads ] ; then echo "number of threads not specified, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping!" ; exit 0 ; fi

## Additional arguments
suffix="" # sample suffix in $dir
export cpath="${ass}/SAMPLE.dipspades/extracted_reads_SAMPLE.fastq.LOCUS.ids.spades/dipspades/consensus_contigs.fasta" # path to raw contigs. SAMPLE and LOCUS will be replaced using regex.
export fpath="${dir}/SAMPLE${suffix}/SAMPLE.LOCUS.bestScore.fasta" # where best contigs / supercontigs are / will be written. SAMPLE and LOCUS will be replaced using regex.

## Export
export ass=$ass
export dir=$dir
export minaln=$minaln
export minnormscore=$minnormscore

## Define function
doCombine()
{
	sample=$1
	echo $sample
	
	combine.contigs.R ${dir}/${sample}${suffix} ${cpath} ${fpath} TRUE TRUE ${minaln} ${minnormscore}
}

## Execute function
export -f doCombine
cat $sfile | parallel -j $threads doCombine

## Finish
echo
echo "All samples processed!"
echo

