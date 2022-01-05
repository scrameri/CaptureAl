#!/bin/bash

#BSUB -J "fastqc[1-192]%30"                         ## run %Y parallel jobs
#BSUB -W 4:00                                       ## 4 hours for each job 
#BSUB -n 1                                          ## 1 core per job -> n * Y cores in parallel
#BSUB -R "rusage[mem=250]"                          ## mem in MB: faster queue if mem < 30000 -> R * Y MB of memory in parallel

## Usage
# bsub < bsub.fastqc.sh

## Needs
# fastqc, java, $sfile (file with sample basename to be assessed, relative to current dir), (optional: adapter.txt)

# load modules
module load gcc/4.8.2 gdc java/1.8.0_73 fastqc/0.11.4

# arguments
in=$(pwd)
out="${in}/fastqc"
sfile="samples.txt"
adapter="/cluster/work/gdc/people/crameris/Dalbergia/uce/seq-qualfiltered/adapters.txt"

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
#if [ ! -f ${adapter} ] ; then echo "ERROR: adapter file <${adapter}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${in} ] ; then echo "ERROR: input directory (reads) <${in}> not found, stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# create output direcory if needed
if [ ! -d ${out} ] ; then mkdir ${out} ; fi

# run job in loop for every read file per sample 
for file in $(ls -1d ${in}/${name}* )
do
	# custom adapter
	#fastqc -t 2 -o $out $file -a $adapter {} 

	# default adapter
	fastqc -t 2 -o ${out} ${file}
done
