#!/bin/bash

#BSUB -J "EXT[1-192]%10"
#BSUB -R "rusage[mem=2500]"          # 2500 should be more than enough for most samples, use 5000 for failed jobs
#BSUB -n 1
#BSUB -W 4:00

## Usage: bsub < bsub.extract.readpairs.sh # submit from ${readdir}

# load modules
module load gdc samtools/1.10 perl/5.16.3

# arguments
sfile="samples.txt"
reg="loci.txt"
Q=10
mapdir=$(pwd) # absolute path needed
readdir="/cluster/scratch/crameris/seq-qualfiltered/$(basename ${mapdir})"
out="/cluster/scratch/crameris/extracted.Q${Q}/$(basename ${readdir})"
bampath="${mapdir}/SAMPLE/SAMPLE.bwa-mem.sorted.Q${Q}.bam" # .nodup.bam

# additional arguments
ext1=".trim1.fastq.gz"
ext2=".trim2.fastq.gz"

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${reg} ] ; then echo "ERROR: locus file <${reg}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${mapdir} ] ; then echo "ERROR: mapping directory <${mapdir}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${readdir} ] ; then echo "ERROR: reads directory <${readdir}> not found, stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# create output directory
if [ ! -d $(dirname ${out}) ] ; then mkdir $(dirname ${out}) ; fi
if [ ! -d ${out} ] ; then mkdir ${out} ; fi
#if [ ! -d ${out}/logs ] ; then mkdir ${out}/logs ; fi

# copy regions and samples
# if reg does not exist, assume all regions in reference *fasta
cp ${reg} ${out}
cp ${sfile} ${out}

# set path to BAM file
bam=$(echo "${bampath}" | sed -e "s/SAMPLE/$name/g")
	
# create subdirectory for each sample
if [ ! -d ${out}/${name} ] ; then mkdir ${out}/${name} ; fi

# Create log file
logfile="${out}/${name}/doExtract.log"

echo "=================================================================" > ${logfile}
echo "========================= doExtract LOG =========================" >> ${logfile}
echo "=================================================================" >> ${logfile}
echo " " >> ${logfile}
echo "starting time:                $(zdump CET)" >> ${logfile}
echo "working directory:            $(pwd)" >> ${logfile}
echo "sample file:                  ${out}/$(basename ${sfile})" >> ${logfile}
echo "region file:                  ${out}/$(basename ${reg})" >> ${logfile}
echo "forward reads used:           ${readdir}/${name}${ext1}" >> ${logfile}
echo "reverse reads used:           ${readdir}/${name}${ext2}" >> ${logfile}
echo "mapped reads used:            ${mapdir}" >> ${logfile}
echo "path to BAM file:             ${bam}" >> ${logfile} 
echo "quality threshold used:       ${Q}" >> ${logfile}
#echo "number of threads used:       ${threads}" >> ${logfile}

# extract mapped reads from .bam file, grep the header (machine-specific grepping!),
# cut the first part of the header, and add a @ at the beginning using a RegExpr
# The @ ensures that the header is identical to the one used in the fastq file
for region in $(cat ${out}/$(basename ${reg}))
   do
   samtools view ${bam} "$region" |cut -f1 |awk '{print "@"$0}' |sort |uniq > ${out}/${name}/${region}.ids
done
			
# create pipes in order to unzip the fastq.gz files on the fly
cd ${out}/${name}
ls -1 *.ids > ${name}.readID

mkfifo ${name}.trim1.fastq
mkfifo ${name}.trim2.fastq
zcat ${readdir}/${name}${ext1} > ${name}.trim1.fastq &
zcat ${readdir}/${name}${ext2} > ${name}.trim2.fastq &

# extract reads using Stefan Zoller's .pl script
extract-reads-from-fastq.pl -f ${name}.trim1.fastq -r ${name}.readID > /dev/null # > ../logs/${name}.extract.reads.log 2> ../logs/${name}.extract.reads.err
extract-reads-from-fastq.pl -f ${name}.trim2.fastq -r ${name}.readID > /dev/null # > ../logs/${name}.extract.reads.log 2> ../logs/${name}.extract.reads.err

# clean up
/bin/rm -f *.ids
/bin/rm -f ${name}.readID
/bin/rm ${name}.trim1.fastq
/bin/rm ${name}.trim2.fastq

# sample finish time
echo "Finish time:                  $(zdump CET)" >> $(basename ${logfile})
echo " " >> $(basename ${logfile})

# tar and compress
cd ${out}
tar -czf ${name}.tar.gz ${name}
/bin/rm -rf ${name}
