#!/bin/bash

## Usage: extract.readpairs.sh -q <qualfiltereddir> -m <mappingdir> -Q <mapping quality threshold> -t <threads>

# -q absolute path to folder with quality-filtered reads
# -m absolute path to folder with mapping dirs
# -t number of threads used

## Needs: samtools, extract-reads-from-fastq.pl

## Authors: Simon Crameri (ETHZ), Stefan Zoller (GDC)

## NOTE: best to start this from a local scratch

## Define arguments
while getopts q:m:Q:t: opts
do
        case "${opts}"
        in
                q) qualfiltereddir=${OPTARG};;
                m) mappingdir=${OPTARG};;
		Q) Q=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done

## Check arguments
if [ ! $qualfiltereddir ] ; then echo "absolute path to quality-filtered reads (-q option) required, stopping." ; exit 0 ; fi
if [ ! $mappingdir ] ; then echo "absolute path to mapped reads (-m option) required, stopping." ; exit 0 ; fi
if [ ! $Q ] ; then echo "mapping quality parameter (-Q option) not set, setting to Q = 10." ; Q=10 ; fi
if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, setting to -t 4." ; threads=4 ; fi

## Create target directory
results=$(basename $mappingdir)
mkdir ${results}
mkdir ${results}/logs

cd ${results}

## Get regions
grep '^>' ${mappingdir}/*.fasta | cut -d '>' -f2 > regions
cp ${mappingdir}/samples.txt .

## Create log file
echo "=================================================================" > doExtract.log
echo "========================= doExtract LOG =========================" >> doExtract.log
echo "=================================================================" >> doExtract.log
echo " " >> doExtract.log
echo "Starting time:                $(zdump MEC)" >> doExtract.log
echo "quality-filtered reads used:  ${qualfiltereddir}" >> doExtract.log
echo "mapped reads used:            ${mappingdir}" >> doExtract.log
echo "quality threshold used: 	    ${Q}" >> doExtract.log 
echo "Number of threads used:       ${threads}" >> doExtract.log


## Define the function
doExtract()
{ 
	sample=$1
	echo "sample: $sample"
	
	
	# create subdirectory for each sample
	mkdir ${sample}.targets
	cd ${sample}.targets
	

	# extract mapped reads from .bam file, grep the header (machine-specific grepping!),
	# cut the first part of the header, and add a @ at the beginning using a RegExpr
	# The @ ensures that the header is identical to the one used in the fastq file
	for region in `cat ../regions`
    	do
		samtools view ${mappingdir}/${sample}/${sample}.bwa-mem.mapped.Q${Q}.sorted.bam "$region" |cut -f1 |awk '{print "@"$0}' |sort |uniq > ${region}.ids
	done
	ls -1 *.ids > ${sample}.readID.files
	
	
	# create pipes in order to unzip the fastq.gz files on the fly
	mkfifo ${sample}.trim1.fastq
	mkfifo ${sample}.trim2.fastq
	zcat ${qualfiltereddir}/${sample}.trim1.fastq.gz > ${sample}.trim1.fastq &
	zcat ${qualfiltereddir}/${sample}.trim2.fastq.gz > ${sample}.trim2.fastq &
	
		
	# extract reads:  needs mem=80000
	extract-reads-from-fastq.pl -f ${sample}.trim1.fastq -r ${sample}.readID.files > ../logs/${sample}.extract.reads.log 2> ../logs/${sample}.extract.reads.err
	extract-reads-from-fastq.pl -f ${sample}.trim2.fastq -r ${sample}.readID.files > ../logs/${sample}.extract.reads.log 2> ../logs/${sample}.extract.reads.err

    cd ..
}

export qualfiltereddir=$qualfiltereddir
export mappingdir=$mappingdir
export Q=$Q
export results=$results
export -f doExtract

cat samples.txt | parallel  -j $threads doExtract


## Sample Finish Time
echo "Finish time:                  $(zdump MEC)" >> doExtract.log
echo " " >> doExtract.log

echo ""
echo "All samples processed."

cd ../ 
