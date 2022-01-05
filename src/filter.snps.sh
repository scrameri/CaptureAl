#!/bin/bash

## Description
# SNP filtering
# based on: John Puritz dDocent pipeline: https://github.com/jpuritz/dDocent/blob/master/tutorials/Filtering%20Tutorial.md

## Usage: filter.snps.sh -v <vcfs.vcf> -r <reference.fasta> -n <output prefix>

## Needs: plot.out.imiss.R if plot.mis='true'

## Load modules (on EULER cluster)
#module load gcc/4.8.2 gdc perl/5.18.4
#module load gcc/4.8.2 gdc vcflib/1.0.1
#module load gcc/4.8.2 gdc python/2.7.11 java/1.8.0_73 perl/5.18.4 freebayes/0.9.20 trimmomatic/0.35 bwa/0.7.12 cd-hit/4.6.4 bedtools/2.25 vcflib/1.0.0 docent/2.12
#module load vcftools/0.1.15
#module load r/3.1.2


## Define arguments
while getopts v:r:n: opts
do
        case "${opts}"
        in
                v) raw=${OPTARG};;
                r) ref=${OPTARG};;
                n) vcfname=${OPTARG};;
    	 esac
done


## Additional arguments
rawname=$(basename $raw .vcf)
refname=$(basename $ref .fasta)
if [ ! $vcfname ] ; then vcfname=$rawname ; fi
plotmis='true'


## Remove individuals if specified in rm.indv
if [ -e rm.indv ] ; then vcftools --vcf $raw --remove rm.indv --recode --recode-INFO-all --out step0 ; else ln -s $raw step0.recode.vcf ; fi


## Get invariant of raw sites
vcftools --vcf $raw --freq --out raw
cat raw.frq | sort -k3,3n -k1,1 -k2,2n | awk '$3 == "1" { print $0 }' > raw.frq.inv
echo "#CHROM POS" > raw.inv
cut -f1,2 raw.frq.inv >> raw.inv
vcftools --vcf step0.recode.vcf --positions raw.inv --recode --recode-INFO-all --out step0.inv


## Get qualitiy of raw sites
egrep -v "^#" $raw | cut -f 6 > raw.qual


## Filter for mapping quality and minimum depth
vcftools --vcf step0.recode.vcf --mac 1 --minDP 3 --minQ 30 --recode --recode-INFO-all --out step1


## Filter individuals that did not sequence well
#--missing-indv Generates a file reporting the missingness on a per-individual basis. The file has the suffix ".imiss".
vcftools --vcf step1.recode.vcf --missing-indv

if [ $plotmis='true' ]
then
	Rscript plot.out.imiss.R                        ##>##>## decide on threshold based on histogram
	mawk '$5 > 1' out.imiss | cut -f1 > lowDP.indv  ##>##>## set threshold here  ## 0 low depth individuals
fi


## Filter for minimum mean depth
#--maf--maf <float> Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value. Allele frequency is defined as the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that site.
#--min-meanDP <float> Includes only sites with mean depth values (over all included individuals) greater than or equal to the "--min-meanDP" value. These options require that the "DP" FORMAT tag is included for each site.
vcftools --vcf step1.recode.vcf --min-meanDP 10 --recode --recode-INFO-all --out step2 


## Decompose complex variants (MNPs) to phased SNPs plus INDELS
vcfallelicprimitives step2.recode.vcf --keep-info --keep-geno > prim.vcf


## Remove indels
#--remove-indels exclude sites that contain an indel. For these options "indel" means any variant that alters the length of the REF allele.
vcftools --vcf prim.vcf --remove-indels --recode --recode-INFO-all --out filtered.SNPs


## Change output name
nsnp=$(mawk '!/#/' filtered.SNPs.recode.vcf | wc -l)
nsnpindels=$(mawk '!/#/' prim.vcf | wc -l)
rawname="$vcfname$nsnp"
mv filtered.SNPs.recode.vcf $rawname.filtered.vcf

rawname2="$vcfnamesnps_indels_$nsnpindels"
mv prim.vcf $rawname2.filtered.vcf


## Rm large intermediate files
rm step*.vcf

