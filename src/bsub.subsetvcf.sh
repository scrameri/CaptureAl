#!/bin/bash

#BSUB -J "SubVCF[1-48]%48"
#BSUB -R "rusage[mem=2500]"
#BSUB -n 1
#BSUB -W 4:00
#BSUB -R "rusage[scratch=20000]"

## Usage: bsub < bsub.subsetvcf.sh

## Description subset a .vcf file by regions and convert it to adegenet::genind and adegenet::DNAbin

## Value: output folder for each locus with one .rda (adegenet objects), one .fasta and one .log file

# load modules
module load samtools/1.2 r/4.0.2

# arguments
vcf="Dalbergia_1197_reg1-171_snps_noprim_341565.filtered.vcf"
ref="/cluster/scratch/crameris/snp.calling/dalbergia.1198.2396loci/consDalbergia_4c_2396.fasta"
batches="batches"

# additional arguments
out="$(basename ${vcf} .vcf)_perlocus"

# check arguments
if [ ! -f ${vcf} ] ; then echo "ERROR: .vcf file <${vcf}> not found, stopping" ; exit 0 ; fi

# define job run variables
IDX=${LSB_JOBINDEX}
regions=${batches}/batch_${IDX}.txt

# create output directory
if [ ! -d ${out} ] ; then mkdir ${out} ; fi

## Subset VCF
# For each parallel job, copy the .vcf to the compute node scratch $TMPDIR
# this will prevent multiple accessions to the file in parallel jobs since the .vcf is 
# copied and accessed separately for each job
vcfbase=$(basename ${vcf} .vcf)
cp ${vcf} ${TMPDIR}/${vcfbase}_${IDX}.vcf

# loop through <n> regions to reduce the amount of times the .vcf needs to be copied
for region in $(cat ${regions})
do
	vcfout="${out}/${region}.vcf"
	
	## grep
	query="^(#|${region}[[:space:]])"
	grep -E ${query} ${TMPDIR}/${vcfbase}_${IDX}.vcf > ${vcfout}
	
	## bcftools (seems a bit slower than grep, and outputs the same variants but 3 additional header lines, and .:.:.:. instead of .\t.\t.)
	#bcftools view ${TMPDIR}/${vcfbase}_${IDX}.vcf --targets ${region} --threads 2 -O v > ${vcfout}
	
	## vcf2gi
	Rscript vcf2gi.R ${vcfout} ${ref} > ${vcfout}.log
	mv ${vcfout}.log ${out}/${region}/vcf2gi.log
	
	## remove .vcf
	#/bin/rm ${out}/${region}/${region}.vcf
done
