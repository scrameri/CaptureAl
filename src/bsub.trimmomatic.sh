#!/bin/bash

#BSUB -J "trimPE[1-192]%10"				# 10 in parallel
#BSUB -R "rusage[mem=3000]"				# 3000 in most cases, 5000 for large fastq files that fail
#BSUB -n 1
#BSUB -W 4:00

## Usage: bsub < bsub.trimmomatic.sh # run from raw reads directory

## Needs: trimmomatic, java, input directory, adapters.txt, (optional samples.txt for a subset of sampels)

# load modules
module load gcc/4.8.2 gdc java/1.8.0_73 trimmomatic/0.35

# arguments
in=$(pwd) # input directory: one line in $sfile & $ext1 or $ext2 must be the path to forward or reverse reads relative to $in
out="/cluster/home/crameris/scratch/seq-qualfiltered/$(basename ${in})" # output directory: will be newly created if inexistent
adapter="/cluster/work/gdc/people/crameris/Dalbergia/uce/seq-qualfiltered/adapters.txt"
sfile="samples.txt"

# create sample list if needed
if [ ! -f ${sfile} ]
then 
	find ${in} -maxdepth 1 -regextype egrep -regex '.*[_.]R1[_.].*' |sed 's!.*/!!' |sed 's/[.][/]//' |sed 's/[_.]R1[_.].*//' |sort |uniq > ${sfile}
fi

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${adapter} ] ; then echo "ERROR: adapter file <${adapter}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${in} ] ; then echo "ERROR: input directory (raw reads) <${in}> not found, stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# create output directory if needed
if [ ! -d $(dirname ${out}) ] ; then mkdir $(dirname ${out}) ; fi
if [ ! -d ${out} ] ; then mkdir ${out} ; fi
if [ ! -d ${out}/logs ] ; then mkdir ${out}/logs ; fi

# copy sample file
cp ${sfile} ${out}

# get read pairs as arrays
reads1=($(ls -1d ${in}/${name}* |sed 's!.*/!!' | grep -E '.*[_.]R1[_.].*'))
reads2=($(ls -1d ${in}/${name}* |sed 's!.*/!!' | grep -E '.*[_.]R2[_.].*'))
len=$(echo ${#reads1[@]})

# run job in loop for every read pair files per sample 
for (( c=0; c<$len; c++ ))
do
	
	# input
	read1=${in}/$(echo ${reads1[c]})
	read2=${in}/$(echo ${reads2[c]})
	
	# output
	trim1=${out}/$(echo $(basename $read1) | sed -E 's/(.*)[_.]R1[_.](.*)/\1.trim1.\2/')
	trim2=${out}/$(echo $(basename $read2) | sed -E 's/(.*)[_.]R2[_.](.*)/\1.trim2.\2/')
	
	unp1=${out}/logs/$(echo $(basename $read1) | sed -E 's/(.*)[_.]R1[_.](.*)/\1.U1.\2/')
	unp2=${out}/logs/$(echo $(basename $read2) | sed -E 's/(.*)[_.]R2[_.](.*)/\1.U2.\2/')
	
	log=${out}/logs/$(echo $(basename $read1) | sed -E 's/(.*)[_.]R1[_.](\w+)(.f.*q.gz)/\1_\2.log/')
	err=${out}/logs/$(echo $(basename $read1) | sed -E 's/(.*)[_.]R1[_.](\w+)(.f.*q.gz)/\1_\2.err/')
	
	# trim default
	#trimmomatic PE -threads 2 -phred33 ${read1} ${read2} ${trim1} ${unp1} ${trim2} ${unp2} ILLUMINACLIP:${adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 > ${log}  2> ${err}

	# trim custom
	trimmomatic PE -threads 2 -phred33 ${read1} ${read2} ${trim1} ${unp1} ${trim2} ${unp2} ILLUMINACLIP:${adapter}:2:20:10:10:true LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50 > ${log} 2> ${err}

done
