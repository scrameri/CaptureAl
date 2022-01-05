#!/bin/bash

#BSUB -J "COV[1-192]%80"
#BSUB -R "rusage[mem=300]"          # 300 is sufficient
#BSUB -n 1
#BSUB -W 4:00

## Usage
# bsub < bsub.get.coverage.stats.sh # from mapping dir

## Needs
# samtools, R, $sfile, sample mapping directories with bwamem output

# load modules
module load gdc samtools/1.2 bedtools/2.28.0 r/3.1.2

# arguments
in=$(pwd)
Q=10
maxcov=100
sfile="samples.txt"

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${in} ] ; then echo "ERROR: input directory (mapping directory) <${in}> not found, stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# set input files
ref=$(echo ${in}/${name}/*.fasta)
bamfileQ="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.bam"
bamfileNODUP="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.cov${maxcov}.bam"
flagstats="${in}/${name}/${name}.flagstats.txt"
flagstatsQ="${in}/${name}/${name}.Q${Q}.flagstats.txt"
flagstatsNODUP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.flagstats.txt"

# set output files
mstatfile="${in}/${name}/${name}.mapstats.txt"
mstatfileQ="${in}/${name}/${name}.Q${Q}.mapstats.txt"
mstatfileNODUP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.mapstats.txt"

# check paths
if [ ! -f ${ref} ] ; then echo "ERROR: reference .fasta file <${ref}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${bamfileQ} ] ; then echo "ERROR: .bam file <${bamfileQ}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${bamfileNODUP} ] ; then echo "ERROR: .bam file <${bamfileNODUP}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${flagstats} ] ; then echo "ERROR: flagstats file <${flagstats}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${flagstatsQ} ] ; then echo "ERROR: flagstatsQ file <${flagstatsQ}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${flagstatsNODUP} ] ; then echo "ERROR: flagstatsNODUP file <${flagstatsNODUP}> not found, stopping" ; exit 0 ; fi

# get mapping stats
# from samtools view: get number of unique read IDs (pairs + singletons)
nidQ=$(samtools view ${bamfileQ} | cut -f1 | sort | uniq | wc -l)
nidNODUP=$(samtools view ${bamfileNODUP} | cut -f1 | sort | uniq | wc -l)
	
# from flagstats: get number of mapped reads (fwd + rev)
tot=$(cat ${flagstats} | grep "in total" | cut -f1 -d' ')
pp=$(cat ${flagstats} | grep "properly paired" | cut -f1 -d' ')
si=$(cat ${flagstats} | grep "singletons" | cut -f1 -d' ')
percpp=$(awk -v a="$tot" -v b="$pp" 'BEGIN {print b/a*100}')

totQ=$(cat ${flagstatsQ} | grep "in total" | cut -f1 -d' ')
ppQ=$(cat ${flagstatsQ} | grep "properly paired" | cut -f1 -d' ')
siQ=$(cat ${flagstatsQ} | grep "singletons" | cut -f1 -d' ')
percppQ=$(awk -v a="$totQ" -v b="$ppQ" 'BEGIN {print b/a*100}')

totNODUP=$(cat ${flagstatsNODUP} | grep "in total" | cut -f1 -d' ')
ppNODUP=$(cat ${flagstatsNODUP} | grep "properly paired" | cut -f1 -d' ')
siNODUP=$(cat ${flagstatsNODUP} | grep "singletons" | cut -f1 -d' ')
percppNODUP=$(awk -v a="$totNODUP" -v b="$ppNODUP" 'BEGIN {print b/a*100}')

echo -e "sample\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${mstatfile}
echo -e "${name}\t${tot}\t${pp}\t${si}\t${percpp}"  >> ${mstatfile}

echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${mstatfileQ}
echo -e "${name}\t${nidQ}\t${totQ}\t${ppQ}\t${siQ}\t${percppQ}"  >> ${mstatfileQ}

echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${mstatfileNODUP}
echo -e "${name}\t${nidNODUP}\t${totNODUP}\t${ppNODUP}\t${siNODUP}\t${percppNODUP}"  >> ${mstatfileNODUP}

# calculate coverage per locus
# write mapped length and average coverage per region to *.coverage.txt
# by default, only considers regions with properly-paired reads (argument 3)
# by default, calculates coverages region by region to save memory in parallel mode (argument 4)	
get.coverage.stats.R ${bamfileQ} ${ref}
get.coverage.stats.R ${bamfileNODUP} ${ref}

# get coverage stats above threshold
covfileQ="${in}/${name}/${name}.Q${Q}.coverage.txt"
covfileNODUP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.coverage.txt"
covtabQ="${in}/${name}/${name}.Q${Q}.covtab.txt"
covtabNODUP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.covtab.txt"
/bin/rm -f ${covtabQ}
/bin/rm -f ${covtabNODUP}

for i in 1 2 3 4 5 6 7 8 9 10 15 20 30 40 50 100 200 300 400 500 600 700 800 900 1000 10000
do 
	echo -n -e  "coverage >= $i\t" >> ${covtabQ}
	cat ${covfileQ} | tail -n +2 | awk '{if($4 >='$i' ) print $0}' |wc -l >> ${covtabQ}

	echo -n -e  "coverage >= $i\t" >> ${covtabNODUP}
	cat ${covfileNODUP} | tail -n +2 | awk '{if($4 >='$i' ) print $0}' |wc -l >> ${covtabNODUP}
done
