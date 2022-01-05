#!/bin/bash

## Usage: trim.fastq.sh -s <sample file> -a <adapter file> -r <directory to raw reads> -x <_R1.fastq.gz,_R2.fastq.gz> -o <output directory> -t <number of threads>

## Arguments
# -s	sample file containing one field (w/o header) with the sample base names
# -a 	adapter file in FASTA format containing adapter sequences
# -r	path to raw reads: will look in this folder for <${sample}${suffix1}> raw reads [DEFAULT: current working directory ]
# -x	suffixes for fwd and rev reads: comma-separated, e.g. '_R1.fastq.gz,_R2.fastq.gz' if the fwd raw read file is ${sample}_R1.fastq.gz and the rev raw reads file is in ${sample}_R2.fastq.gz
# -o	path to output directory [DEFAULT: current working directory ]
# -t 	threads (should not exceed 30)

## Needs: 
# java, trimmomatic

## Define arguments
while getopts s:a:r:x:o:t: opts
do
        case "${opts}"
		in
				s) sfile=${OPTARG};;
				a) adapterfile=${OPTARG};;
				r) readpath=${OPTARG};;
				x) suffix=${OPTARG};;
				o) odir=${OPTARG};;
				t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file not specified (-s option), stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! $adapterfile ] ; then echo "adapter file not specified (-a option), stopping" ; exit 0 ; fi
if [ ! -f $adapterfile ] ; then echo "adapter file <$adapter> not found, stopping" ; exit 0 ; fi
if [ ! $readpath ] ; then echo "path to raw reads (-r option) not provided, using current directory" ; readpath=$(pwd) ; fi
if [ ! $odir ] ; then echo "path to output directory (-o option) not provided, using current directory" ; odir=$(pwd) ; fi
if [ ! -d $readpath ] ; then echo "path to raw reads <$readpath> not found, stopping" ; exit 0 ; fi
if [ ! $suffix ] ; then echo "sample suffixes not found, assuming '_R1.fastq.gz,_R2.fastq.gz'" ; suffix="_R1.fastq.gz,_R2.fastq.gz" ; fi
if [ ! $threads ] ; then echo "number of threads not specified, using threads=16" ; threads=16 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads (-t option) must not exceed 30, stopping." ; exit 0 ; fi

## Internal arguments
export suffix1=$(echo ${suffix} | cut -d "," -f1)
export suffix2=$(echo ${suffix} | cut -d "," -f2)

## Create output directory
if [[ ! -d $odir ]] ; then mkdir $odir ; fi

## Create LOG file
echo "=================================================================" > ${odir}/trimming.log
echo "========================== Trimming LOG =========================" >> ${odir}/trimming.log
echo "=================================================================" >> ${odir}/trimming.log
echo " " >> ${odir}/trimming.log
echo "Starting time:                $(zdump MEC)" >> ${odir}/trimming.log
echo "sample file:                  ${sfile}" >> ${odir}/trimming.log
echo "adapter file:                 ${adapterfile}" >> ${odir}/trimming.log
echo "path to raw reads:            ${readpath}" >> ${odir}/trimming.log
echo "fwd reads suffix:             ${suffix1}" >> ${odir}/trimming.log
echo "rev reads suffix:             ${suffix2}" >> ${odir}/trimming.log
echo "number of threads:            ${threads}" >> ${odir}/trimming.log

## Set JAVA options
ulimit -s unlimited
export _JAVA_OPTIONS="-XX:ParallelGCThreads=2"

## Define trimming function using Trimmomatic
doTrimming() {
	base=$1 # sample basename
	
	read1=${readpath}/${base}${suffix1} # path to fwd reads
	read2=${readpath}/${base}${suffix2} # path to rev reads
	echo ${base} $(basename ${read1}) $(basename ${read2})
	
	# java -Xms1G -Xmx3G -Xss512k -jar /usr/local/trimmomatic/trimmomatic-0.32.jar PE -threads 1
	trimmomatic PE -threads 1 -phred33  $read1  $read2  ${odir}/${base}.trim1.fastq.gz ${odir}/${base}.trim1.unp.fastq.gz ${odir}/${base}.trim2.fastq.gz ${odir}/${base}.trim2.unp.fastq.gz ILLUMINACLIP:${adapterfile}:2:20:10:10:true LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50  >  ${odir}/${base}.trimmo.log 2> ${odir}/${base}.trimmo.err 
}

## Execute function
export readpath=$readpath
export adapterfile=$adapterfile
export odir=$odir
export -f doTrimming 

cat ${sfile} | parallel  -j $threads doTrimming


## Finish
echo "Finish time:                  $(zdump MEC)" >> ${odir}/trimming.log
echo " " >> ${odir}/trimming.log
echo
echo "All samples processed!"
echo
