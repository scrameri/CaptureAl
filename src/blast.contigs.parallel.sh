#!/bin/bash

## Usage: blast.contigs.parallel.sh -s <sample file> -r <reference.fasta> -d <contigs directory> -c <contig fasta name> -e <optional: e-value [1e-04]> -t <optional: number of threads [4]>

## Needs: BLAST+, blast.fasta.seqs.sh

## Define arguments
while getopts s:r:d:c:e:t: opts
do
        case "${opts}"
        in
		s) sfile=${OPTARG};;
		r) ref=${OPTARG};;
		d) blastdir=${OPTARG};;
		c) contigs=${OPTARG};;
		e) evalue=${OPTARG};;
		t) threads=${OPTARG};;

    	esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file not specified (-s option), stopping" ; exit 0 ; fi
if [ ! $ref ] ; then echo "reference fasta file not specified (-r option), stopping" ; exit 0 ; fi
if [ ! $blastdir ] ; then echo "directory with sample subdirs with contigs to blast not specified (-d option), setting to 'blast.contigs'" ; blastdir=blast.contigs ; fi
if [ ! $contigs ] ; then echo "name of contigs fasta in sample subdirs not specified (-c option), stopping" ; exit 0 ; fi

if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! -f $ref ] ; then echo "reference fasta <$ref> not found, stopping" ; exit 0 ; fi
if [ ! -d $blastdir ] ; then echo "directory <$blastdir> not found, stopping" ; exit 0 ; fi

if [ ! $evalue ] ; then echo "evalue not specified, setting to 1e-04" ; evalue=1e-04 ; fi
if [ ! $threads ] ; then echo "number of threads not specified, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments


## Export arguments
export ref=$ref
export blastdir=$blastdir
export evalue=$evalue
export contigs=$contigs

export dir=$(pwd)
export refbase=$(basename $ref .fasta)
export contigbase=$(basename ${contigs} .fasta)

## Be verbose
echo ""
echo "${ref} used as BLAST queries"
echo "${contigs} used as BLAST subjects"
echo "evalue set to ${evalue}"
echo "${threads} threads used"
echo ""

## Define function
doBlast()
{
	sample=$1
	echo $sample
		
	cd ${blastdir}/${sample}
	if [ -f $ref ] ; then unlink $ref ; fi
	ln -s ${dir}/$ref .

	# check that contigs are available
	if [ ! -f ${contigs} ] 
	then
		echo "$contigs not found in ${blastdir}/${sample}, stopping" >> ${dir}/blast.err
		exit 0
	fi

	# create BLAST database of subject sequences $contigs)
	if [ ! -f ${contigs}.nog ]
	then
		makeblastdb -in ${contigs} -parse_seqids -dbtype nucl >/dev/null
	fi

	if [ ! -f ${contigs}.nog ] ; then echo "could not create blast database for $sample, stopping" >> ${dir}/blast.err ; fi
	
	# BLAST using 1 thread
	blast.fasta.seqs.sh -d ${contigs} -q ${ref} -e ${evalue} -t 1 >/dev/null
	
	# rename output
	mv ${refbase}.on.${contigbase}.blast ${refbase}.on.${sample}.${contigbase}.${evalue}.blast
	
	cd ${dir}
}

## Execute function
export -f doBlast
cat ${sfile} | parallel -j $threads doBlast

## Finish
echo ""
echo "All samples processed."
echo ""

