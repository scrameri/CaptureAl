#!/bin/bash

#BSUB -J "BWA[1-192]%30"
#BSUB -R "rusage[mem=4333]"          # 4333 requests 2 x 4333 if -n 2 ; use 13000 % 5 for failed jobs  
#BSUB -n 2
#BSUB -W 4:00
#BSUB -R "rusage[scratch=20000]"

## Usage
# bsub < bsub.bwamem.sh # run from trimmed reads directory

## Needs
# $fas, $sfile, bwa, sambamba, java, python, picard-tools, jvarkit

# load modules
module load gdc gcc/4.8.2 bwa/0.7.17 sambamba/0.8.0 samtools/1.2 java/1.8.0_101 python/3.6.1 java/1.8.0_73 picard-tools/2.23.8 r/3.1.2

# arguments
in=$(pwd) # input directory: one line in $sfile & $ext1 or $ext2 must be the path to forward or reverse reads relative to $in
sortsamrefname="/cluster/project/gdc/shared/tools/jvarkit/dist/sortsamrefname.jar" # jvarkit tool required for capping
biostar154220="/cluster/project/gdc/shared/tools/jvarkit/dist/biostar154220.jar" # jvarkit tool required for capping
out="/cluster/home/crameris/scratch/mapping-reads-to-2396/$(basename ${in})" # output directory: will be newly created if inexistent
fas="/cluster/home/crameris/home/Dalbergia/uce/references/consDalbergia_4c_2396.fasta" # reference fasta file
sfile="samples.txt" # one line in $sfile & $ext1 or $ext2 must be the path to forward or reverse reads

# additional arguments
mode="paired" # paired or unpaired
ext1=".trim1.fastq.gz" # one line in $sfile & $ext1 must be the path to forward reads
ext2=".trim2.fastq.gz" # one line in $sfile & $ext2 must be the path to reverse reads

cpu=4 # processors per job (default: 4)
T=10 # bwa mem: dont output with score lower than T (default 30)
Q=10 # mapping quality (default: 20)

maxcov=500
maxmem='5g'

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${fas} ] ; then echo "ERROR: reference .fasta file <${fas}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${in} ] ; then echo "ERROR: input directory (trimmed reads) <${in}> not found, stopping" ; exit 0 ; fi

if [ ! -f ${sortsamrefname} ] ; then echo "ERROR: jvarkit executable <${sortsamrefname}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${biostar154220} ] ; then echo "ERROR: jvarkit executable <${biostar154220}> not found, stopping" ; exit 0 ; fi

if [ "$T" -gt "$Q" ] ; then echo "high-quality parameter (Q) must be equal or greater than mapping quality parameter (T), stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# create output directories if needed
if [ ! -d $(dirname $out) ] ; then mkdir $(dirname $out) ; fi
if [ ! -d ${out} ] ; then mkdir ${out} ; fi

# create sample subdirectories
if [ ! -d ${out}/${name} ] ; then mkdir ${out}/${name} ; fi

# copy reference.fasta and fix headers if needed (replace ' ' and '|' with '_')
ref="${out}/${name}/$(basename ${fas})"
if [ ! -f ${ref} ] ; then cp ${fas} ${out}/${name} ; fi
fix.fasta.headers.R ${ref} ' |\\|' '_' FALSE

# copy samples file and create regions file
cp ${sfile} ${out}
grep '^>' ${ref} |cut -f2 -d'>' > ${out}/loci.txt

# index reference fasta
if [ ! -f ${ref}.bwt ] ; then bwa index ${ref} ; fi

# create log file
logfile=${out}/${name}/bwa.Q${Q}.log

echo "=================================================================" > ${logfile}
echo "========================== BWA MEM LOG ==========================" >> ${logfile}
echo "=================================================================" >> ${logfile}
echo " " >> ${logfile}
echo "Starting time:                $(zdump CET)" >> ${logfile}
echo "Working directory:            $(pwd)" >> ${logfile}
echo "Input directory:              ${in}" >> ${logfile}
echo "Output directory:             ${out}" >> ${logfile}
echo "Sample file:                  ${sfile}" >> ${logfile}
echo "Reference used for BWA MEM:   ${fas}" >> ${logfile}
#echo "File extension:               ${ext}" >> ${logfile}
echo "Mode:                         ${mode}" >> ${logfile}
echo "Path to sample:               ${in}/${name}${ext1}" >> ${logfile}
if [ $mode = 'paired' ] ; then echo "Path to sample [pair]:        ${in}/${name}${ext2}" >> ${logfile} ; fi
echo "bwa mem quality parameter -T: ${T}" >> ${logfile}
echo "mapping_quality parameter Q:  ${Q}" >> ${logfile}
echo "Processors per job:           ${cpu}" >> ${logfile}
echo "Maximum coverage (capping):   ${maxcov}" >> ${logfile}
echo "Maximum memory (capping)      ${maxmem}" >> ${logfile}

## MAPPING: BWA MEM
# some bwa mem options
# t: threads
# M: flag shorter split hits as secondary (important for picard tools)
# w: band width, essentially gaps longer than this are not found (default 100)
# T: dont output with score lower than T (default 30)
# B: mismatch penalty (default 4)
# E: gap extension penalty (default 1)

# map reads using bwa mem
samfile="${TMPDIR}/${name}.bwa-mem.Q${Q}.sam"

if [ $mode = 'paired' ]
then

	bwa mem -t ${cpu} -M -T ${T} -R "@RG\tID:$name\tSM:$name" ${ref} ${in}/${name}${ext1} ${in}/${name}${ext2} > ${samfile} #2> ${out}/${name}/${name}.bwa-mem.Q${Q}.log

elif [ $mode = 'unpaired' ]
then
	
	bwa mem -t ${cpu} -M -T ${T} -R "@RG\tID:$name\tSM:$name" ${ref} ${in}/${name}${ext1} > ${samfile} #2> ${out}/${name}/${name}.bwa-mem.Q${Q}.log
		
fi

# convert .sam to sorted.bam
bamfile="${TMPDIR}/${name}.bwa-mem.sorted.Q${Q}.bam"
sambamba view -t ${cpu} -S ${samfile} -f bam -o /dev/stdout |sambamba sort /dev/stdin -o /dev/stdout -t ${cpu} -l 0 -m 6GB --tmpdir ${TMPDIR} > ${bamfile}

# mapping flagstats
flagstats="${out}/${name}/${name}.flagstats.txt"
sambamba flagstat ${bamfile} > ${flagstats}

# remove reads with low mapping quality
bamfileQ="${out}/${name}/${name}.bwa-mem.sorted.Q${Q}.bam"
sambamba view -F "mapping_quality >= ${Q}" ${bamfile} -o ${bamfileQ} -t ${cpu} -f bam 

# remove PCR duplicates and write to sample directory
dupstats="${out}/${name}/${name}.dupstats.txt"
bamfileNODUP="${TMPDIR}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam"
picard MarkDuplicates -I ${bamfileQ} -O ${bamfileNODUP} -M ${dupstats} -TMP_DIR ${TMPDIR} -VALIDATION_STRINGENCY LENIENT -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 512 -REMOVE_DUPLICATES true

# cap nodup.bam file
bamfileCAP="${out}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.cov${maxcov}.bam"
java -Xmx${maxmem} -jar ${sortsamrefname} ${bamfileNODUP} |\
java -Xmx${maxmem} -jar ${biostar154220} -n ${maxcov} 2> /dev/null |\
samtools sort -T ${name} -o ${bamfileCAP}

# create index files
sambamba index ${bamfileCAP}

# final mapping flagstats
flagstatsQ="${out}/${name}/${name}.Q${Q}.flagstats.txt"
flagstatsCAP="${out}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.flagstats.txt"
sambamba flagstat ${bamfileQ} > ${flagstatsQ}
sambamba flagstat ${bamfileCAP} > ${flagstatsCAP}

# log finish Time
echo "Finish time:                  $(zdump CET)" >> ${logfile}
echo " " >> ${logfile}

# tar.gz sample subdirectory with mapping results
#tar -czf ${out}/${name}.tar.gz ${name}
#/bin/rm -rf ${out}/${name}
