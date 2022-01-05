#!/bin/bash

#BSUB -J "TAR[1-192]%30"
#BSUB -R "rusage[mem=500]"         # 500
#BSUB -n 1
#BSUB -W 4:00

# Usage: bsub < bsub.tar.samples.sh # from any directory with untarred / unzipped results (sample directories)

# arguments
sfile="samples.txt"

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# tar and compress
tar -czf ${name}.tar.gz ${name}

# remove original
if [ -d ${name} ] ; then /bin/rm -rf ${name} ; fi
