#!/bin/bash

## Description a .vcf file by regions provided in a region file

## Usage: subset.vcf.by.region.sh -v <.vcf> -r <regions.txt with desired region name(s) [without '>']> -t <number of threads>

## Value: output folder with one .vcf file per region

## Define arguments
while getopts v:r:t: opts
do
        case "${opts}"
        in
                v) vcf=${OPTARG};;
                r) reg=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done

## Check arguments
if [ ! $vcf ] ; then echo "vcf file (-v option) not provided, stopping!" ; exit 0 ; fi
if [ ! -f $vcf ] ; then echo "vcf file <$vcf> does not exist, stopping!" ; exit 0 ; fi
if [ ! $reg ] ; then echo "region file (-r option) not provided, stopping!" ; exit 0 ; fi
if [ ! -f $reg ] ; then echo "region file <$reg> does not exist, stopping!" ; exit 0 ; fi
if [ ! $threads ] ; then echo "number of threads not provided, setting to 4" ; threads=4 ; fi

## Export
export vcf=$vcf
export vcfbase=$(basename $vcf .vcf)
if [ -d $vcfbase ] ; then echo "output directory <$vcfbase> already exists, stopping!" ; exit 0 ; else mkdir $vcfbase ; fi

## Define function
doSubsetVcf()
{
	region=$1
	echo $region

	vcfout="${region}.vcf"
	query="^(#|${region}[[:space:]])"
    grep -E ${query} $vcf > ${vcfbase}/$vcfout

}

## Execute function
export -f doSubsetVcf
cat $reg | parallel -j $threads doSubsetVcf

## Finish
echo ""
echo "All regions processed!"
echo ""

