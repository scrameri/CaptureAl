#!/bin/bash
#BSUB -J "mergeBAM"
#BSUB -R "rusage[mem=4000]"
#BSUB -n 4
#BSUB -W 24:00

###Nik Zemp, GDC 15.03.2021, version 0.1
####work on the scratch and keep the temporary bam files here. You can always merge them again.

# arguments
threads=8
bamlist="bamlist.txt"
ofile="merged.bam"

# load modules
module load gdc gcc/4.8.2 perl/5.18.4 samtools/1.10

# merge bamfiles
samtools merge -b ${bamlist} --threads ${threads} ${ofile}

# infex output .bam
samtools index ${ofile}
