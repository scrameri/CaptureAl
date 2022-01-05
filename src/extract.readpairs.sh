#!/bin/bash

# best to start this from a local scratch

## Usage: extract.readpairs.sh -s <samples> -l <regions to extract reads from> -d <qualfiltereddir> -m <mapdir> -o <output directory> -Q <mapping quality threshold> -b <path to BAM file> -t <threads>

## Needs: extract-reads-from-fastq.pl

# -s sample file
# -l locus file
# -d absolute path to folder with quality-filtered reads
# -m absolute path to folder with mapping dirs
# -o output directory (created if inexistent)
# -Q mapping quality
# -b path to BAM file (<SAMPLE> can be part of the string and will be replaced by the actual sample using regex)
# -t number of threads used

## Define arguments
while getopts s:l:d:m:o:Q:b:t: opts
do
        case "${opts}"
        in
       		s) sfile=${OPTARG};;
       		l) reg=${OPTARG};;
        	d) readdir=${OPTARG};;
        	m) mapdir=${OPTARG};;
        	o) odir=${OPTARG};;
			Q) Q=${OPTARG};;
			b) bampath=${OPTARG};;
			t) threads=${OPTARG};;
    	 esac
done

## Check arguments
if [ ! ${sfile} ] ; then echo "sample file (-s option) not specified, assuming <samples.txt>." ; sfile="samples.txt" ; fi
if [ ! ${reg} ] ; then echo "region file (-l option) required, stopping." ; exit 0 ; fi
if [ ! -f ${sfile} ] ; then echo "sample file <${sfile}> not found, stopping." ; exit 0 ; fi
if [ ! -f ${reg} ] ; then echo "region file <${reg}> not found, stopping." ; exit 0 ; fi

if [ ! ${readdir} ] ; then echo "absolute path to quality-filtered reads (-d option) required, stopping." ; exit 0 ; fi
if [ ! ${mapdir} ] ; then echo "absolute path to mapped reads (-m option) required, stopping." ; exit 0 ; fi
if [ ! $odir ] ; then echo "output directory (-o option) not specified, using <seq-extracted>" ; odir="seq-extracted" ; fi
if [ ! -d ${readdir} ] ; then echo "folder with quality-filtered reads <${readdir}> not found, stopping." ; exit 0 ; fi
if [ ! -d ${mapdir} ] ; then echo "folder with mapped reads <${mapdir}> not found, stopping." ; exit 0 ; fi

if [ ! ${Q} ] ; then echo "mapping quality parameter (-Q option) not set, assuming -Q 10." ; Q=10 ; fi

if [ ! ${bampath} ] ; then echo "path to BAM file not provided (-b option), setting to <${mapdir}/SAMPLE/SAMPLE.bwa-mem.sorted.Q${Q}.nodup.bam>." ; bampath="${mapdir}/SAMPLE/SAMPLE.bwa-mem.sorted.Q${Q}.nodup.bam" ; fi

if [ ! ${threads} ] ; then echo "number of threads (-t option) not specified, using -t 4." ; threads=4 ; fi
if [ "$threads" -gt 30 ] ; then echo "number of threads (-t option) must be between 1 and 30, setting -t 2." ; threads=2 ; fi

## Additional arguments
ext1=".trim1.fastq.gz"
ext2=".trim2.fastq.gz"

## Create target directory
if [ -d ${odir} ] ; then echo "output directory <${odir}> already exists, moving to <${odir}.bak>" ; mv ${odir} ${odir}.bak ; fi
mkdir ${odir}
#mkdir ${odir}/logs

## Copy regions and samples
cp ${reg} ${odir}/loci.txt
cp ${sfile} ${odir}/samples.txt

## Create log file
logfile="${odir}/doExtract.log"

echo "=================================================================" > ${logfile}
echo "========================= doExtract LOG =========================" >> ${logfile}
echo "=================================================================" >> ${logfile}
echo " " >> ${logfile}
echo "starting time:                $(zdump CET)" >> ${logfile}
echo "working directory:            $(pwd)" >> ${logfile}
echo "sample file:                  ${sfile}" >> ${logfile}
echo "region file:                  ${reg}" >> ${logfile}
echo "quality-filtered reads used:  ${readdir}" >> ${logfile}
echo "mapped reads used:            ${mapdir}" >> ${logfile}
echo "path to BAM file:             ${bampath}" >> ${logfile} 
echo "quality threshold used:       ${Q}" >> ${logfile}
echo "number of threads used:       ${threads}" >> ${logfile}

## Define the function
doExtract()
{ 
	# module load gcc/4.8.2 gdc samtools # only used on euler 
	name=$1
	echo $name
	
	# set path to BAM file
	bam=$(echo "${bampath}" | sed -e "s/SAMPLE/$name/g")
	
	# create subdirectory for each sample
	mkdir ${odir}/${name}

	# extract mapped reads from .bam file, grep the header (machine-specific grepping!),
	# cut the first part of the header, and add a @ at the beginning using a RegExpr
	# The @ ensures that the header is identical to the one used in the fastq file
	for region in $(cat ${odir}/loci.txt)
    do
    	samtools view ${bam} "$region" |cut -f1 |awk '{print "@"$0}' |sort |uniq > ${odir}/${name}/${region}.ids

	done
			
	# create pipes in order to unzip the fastq.gz files on the fly
	cd ${odir}/${name}
	ls -1 *.ids > ${name}.readID.files

	mkfifo ${name}.trim1.fastq
	mkfifo ${name}.trim2.fastq
	zcat ${readdir}/${name}${ext1} > ${name}.trim1.fastq &
	zcat ${readdir}/${name}${ext2} > ${name}.trim2.fastq &
	
	# extract reads:  needs mem=80000	
	extract-reads-from-fastq.pl -f ${name}.trim1.fastq -r ${name}.readID.files > /dev/null # > ../logs/${name}.extract.reads.log 2> ../logs/${name}.extract.reads.err
	extract-reads-from-fastq.pl -f ${name}.trim2.fastq -r ${name}.readID.files > /dev/null # > ../logs/${name}.extract.reads.log 2> ../logs/${name}.extract.reads.err

	# clean up
	/bin/rm -f *.ids
	/bin/rm -f ${name}.readID.files
	/bin/rm ${name}.trim1.fastq
	/bin/rm ${name}.trim2.fastq

    cd ../../
}

export odir=${odir}
export readdir=${readdir}
export mapdir=${mapdir}
export Q=${Q}
export bampath=${bampath}
export ext1=${ext1}
export ext2=${ext2}
export -f doExtract

cat ${odir}/samples.txt | parallel  -j $threads doExtract


## Sample Finish Time
echo "Finish time:                  $(zdump CET)" >> ${logfile}
echo " " >> ${logfile}

echo ""
echo "All samples processed."
