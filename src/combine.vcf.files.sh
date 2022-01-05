#!/bin/bash

## Usage: combine.vcf.files.sh <folder with vcfs to be combined> 

# load modules
module load gdc vcflib/1.0.1

# arguments
folder=$1
fname=$(basename $folder)

# additional arguments
#vcffirstheader="~/bin/freebayes/vcflib/scripts/vcffirstheader"
#vcfstreamsort="~/bin/freebayes/vcflib/bin/vcfstreamsort"

cat ${folder}/*.vcf | vcffirstheader | vcfstreamsort -w 1000 | vcfuniq > ${fname}.vcf
