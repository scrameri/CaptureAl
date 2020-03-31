#!/bin/bash

# best to start this from a local scratch

## Usage: extract.readpairs.sh -s <samples> -r <regions to extract reads from> -d <qualfiltereddir> -m <mappingdir> -Q <mapping quality threshold> -t <threads>

## Needs: extract-reads-from-fastq.pl

# -s sample file
# -r locus file
# -d absolute path to folder with quality-filtered reads
# -m absolute path to folder with mapping dirs
# -Q mapping quality
# -t number of threads used

## Define arguments
while getopts s:r:d:m:Q:t: opts
do
        case "${opts}"
        in
        		s) sfile=${OPTARG};;
        		r) reg=${OPTARG};;
                d) qualfiltereddir=${OPTARG};;
                m) mappingdir=${OPTARG};;
				Q) Q=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) required, stopping." ; exit 0 ; fi
if [ ! $reg ] ; then echo "region file (-r option) required, stopping." ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping." ; exit 0 ; fi
if [ ! -f $reg ] ; then echo "region file <$reg> not found, stopping." ; exit 0 ; fi

if [ ! $qualfiltereddir ] ; then echo "absolute path to quality-filtered reads (-q option) required, stopping." ; exit 0 ; fi
if [ ! $mappingdir ] ; then echo "absolute path to mapped reads (-m option) required, stopping." ; exit 0 ; fi
if [ ! -d $qualfiltereddir ] ; then echo "folder with quality-filtered reads <$qualfiltereddir> not found, stopping." ; exit 0 ; fi
if [ ! -d $mappingdir ] ; then echo "folder with mapped reads <$mappingdir> not found, stopping." ; exit 0 ; fi

if [ ! $Q ] ; then echo "mapping quality parameter (-Q option) not set, assuming -Q 10." ; Q=10 ; fi
if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, using -t 4." ; threads=4 ; fi
if [ "$threads" -gt 30 ] ; then echo "number of threads (-t option) must be between 1 and 30, setting -t 4." ; threads=4 ; fi

## Create target directory
results=$(basename $mappingdir)
mkdir ${results}
mkdir ${results}/logs

## Copy regions and samples
cp ${reg} ${results}/regions
cp ${sfile} ${results}/samples.txt
cd ${results}

## Create log file
echo "=================================================================" > doExtract.log
echo "========================= doExtract LOG =========================" >> doExtract.log
echo "=================================================================" >> doExtract.log
echo " " >> doExtract.log
echo "Starting time:                $(zdump MEC)" >> doExtract.log
echo "sample file:                  $sfile" >> doExtract.log
echo "region file:                  $reg" >> doExtract.log
echo "quality-filtered reads used:  ${qualfiltereddir}" >> doExtract.log
echo "mapped reads used:            ${mappingdir}" >> doExtract.log
echo "quality threshold used: 	    ${Q}" >> doExtract.log 
echo "Number of threads used:       ${threads}" >> doExtract.log


## Define the function
doExtract()
{ 
	# module load gcc/4.8.2 gdc samtools # only used on euler 
	sample=$1
	echo $sample
	
	
	# create subdirectory for each sample
	mkdir ${sample}.targets
	cd ${sample}.targets
	

	# extract mapped reads from .bam file, grep the header (machine-specific grepping!),
	# cut the first part of the header, and add a @ at the beginning using a RegExpr
	# The @ ensures that the header is identical to the one used in the fastq file
	for region in $(cat ../regions)
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
