#!/bin/bash

## Usage: extract.fasta.from.tar.parallel.sh -d <directory with .tar directories> -s <sample file with .tar file basenames> -l <locus file with .fasta basenames to extract> -t <threads>

## Define arguments
while getopts d:s:l:t: opts
do
        case "${opts}"
		in
				d) indir=${OPTARG};;	
				s) sfile=${OPTARG};;
				l) lfile=${OPTARG};;				
				t) threads=${OPTARG};;
		esac
done

## Check parameters and set defaults
if [ ! $indir ] ; then echo "directory with .tar directories not provided (-d option), stopping" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "directory <$indir> does not exist, stopping" ; exit 0 ; fi
if [ ! $lfile ] ; then echo "locus file not provided (-l option), stopping" ; exit 0 ; fi
if [ ! -f $lfile ] ; then echo "locus file <$lfile> does not exist, stopping" ; exit 0 ; fi
if [ ! $sfile ] ; then echo "samples file not provided (-s option), stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "samples file <$sfile> does not exist, stopping" ; exit 0 ; fi
if [ ! $threads ] ; then echo threads=1 ; "using 1 threads" ; else echo "using $threads threads" ; fi

## Export arguments
export indir=$indir
export sfile=$sfile
export lfile=$lfile
export suffix=".bestScore.fasta" # .fasta suffix after name provided in $lfile

## Define extractor function
doExtractFromTar()
{
	s=$1
	echo "extracting sample $s"
	
	for l in $(cat $lfile)
	do	
		# echo "  processing locus $l"
		cd $indir
		tar -xf ${s}.tar ${s}/${s}.${l}.bestScore.fasta
		cd ../
	done	
}

## Execute function
export -f doExtractFromTar
cat $sfile | parallel -j $threads doExtractFromTar

## Finish
echo ""
echo "All samples processed."
echo ""
