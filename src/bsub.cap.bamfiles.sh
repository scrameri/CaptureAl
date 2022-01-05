#!/bin/bash

#BSUB -J "CAP[1-192]%15"
#BSUB -R "rusage[mem=4333]"         # 4333 requests 2 x 4333 if -n 2 ; use 13000 % 5 for failed jobs
#BSUB -n 1
#BSUB -W 4:00
#BSUB -R "rusage[scratch=20000]"

## Usage: bsub < bsub.cap.bamfiles.sh # from mapping directory
# This script is also integrated in bsub.bwamem.sh, but can be used to cap .bam files with other ${maxcov} values

# load modules
module load gdc gcc/4.8.2 sambamba/0.8.0 samtools/1.2 bedtools/2.28.0 java/1.8.0_101 java/1.8.0_73 picard-tools/2.23.8 r/3.1.2

# arguments
in=$(pwd)
Q=10
maxcov=50
maxmem='5g'
sfile="samples.txt"
sortsamrefname="/cluster/project/gdc/shared/tools/jvarkit/dist/sortsamrefname.jar"
biostar154220="/cluster/project/gdc/shared/tools/jvarkit/dist/biostar154220.jar"

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${sortsamrefname} ] ; then echo "ERROR: jvarkit exectuable <${sortsamrefname}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${biostar154220} ] ; then echo "ERROR: jvarkit exectuable <${biostar154220}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${in} ] ; then echo "ERROR: input directory (mapping directory) <${in}> not found, stopping" ; exit 0 ; fi
if [ ! ${maxcov} ] ; then echo "MESSAGE: maximum coverage not defined, setting to 1000" ; maxcov=1000 ; fi
if [ ${maxcov} -lt 1 ] ; then echo "ERROR: maximum coverage must be greater than zero, stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# define input file
bamfileQ="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.bam"

if [ ! -f ${bamfileQ} ] ; then echo "input .bam file <${bamfileQ}> not found, stopping" ; exit 0 ; fi

# remove PCR duplicates and write to sample directory
dupstats="${in}/${name}/${name}.dupstats.txt"
bamfileNODUP="${TMPDIR}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam"
picard MarkDuplicates -I ${bamfileQ} -O ${bamfileNODUP} -M ${dupstats} -TMP_DIR ${TMPDIR} -VALIDATION_STRINGENCY LENIENT -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 512 -REMOVE_DUPLICATES true

# cap nodup.bam file
bamfileCAP="${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.cov${maxcov}.bam"
java -Xmx${maxmem} -jar ${sortsamrefname} ${bamfileNODUP} |\
java -Xmx${maxmem} -jar ${biostar154220} -n ${maxcov} 2> /dev/null |\
samtools sort -T ${name} -o ${bamfileCAP}

# index ${bamfileCAP}
samtools index ${bamfileCAP}

# coverage analysis
# set input / output files
flagstatsCAP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.flagstats.txt"
mstatCAP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.mapstats.txt"
ref=$(echo ${in}/${name}/*.fasta)

# get mapping stats
sambamba flagstat ${bamfileCAP} > ${flagstatsCAP}

# from samtools view: get number of unique read IDs (pairs + singletons)
nidCAP=$(samtools view ${bamfileCAP} | cut -f1 | sort | uniq | wc -l)

# from flagstats: get number of mapped reads (fwd + rev)
totCAP=$(cat ${flagstatsCAP} | grep "in total" | cut -f1 -d' ')
ppCAP=$(cat ${flagstatsCAP} | grep "properly paired" | cut -f1 -d' ')
siCAP=$(cat ${flagstatsCAP} | grep "singletons" | cut -f1 -d' ')
percppCAP=$(awk -v a="$totCAP" -v b="$ppCAP" 'BEGIN {print b/a*100}')

echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${mstatCAP}
echo -e "${name}\t${nidCAP}\t${totCAP}\t${ppCAP}\t${siCAP}\t${percppCAP}"  >> ${mstatCAP}

# calculate coverage per locus
get.coverage.stats.R ${bamfileCAP} ${ref}

# get coverage stats above threshold
covfileCAP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.coverage.txt"
covtabCAP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.covtab.txt"
/bin/rm -f ${covtabCAP}

for i in 1 2 3 4 5 6 7 8 9 10 15 20 30 40 50 100 200 300 400 500 600 700 800 900 1000 10000
do 
	echo -n -e  "coverage >= $i\t" >> ${covtabCAP}
	cat ${covfileCAP} | tail -n +2 | awk '{if($4 >='$i' ) print $0}' |wc -l >> ${covtabCAP}
done
