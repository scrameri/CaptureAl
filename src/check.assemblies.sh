#!/bin/bash

## Usage: check.assemblies.sh -s <samples.txt> -l <locus file> -d <input directory> -f <suffix of sample subdirectory> -x <suffix of assembly subfolder in sample subdirectory>

## Get arguments
while getopts s:l:d:f:x: opts
do
        case "${opts}"
        in
        	s) sfile=${OPTARG};;
		l) lfile=${OPTARG};;
        	d) indir=${OPTARG};;
		f) samsuf=${OPTARG};;
		x) suffix=${OPTARG};;
        esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file not specified (-s option), stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! $lfile ] ; then echo "locus file not specified (-l option), stopping" ; exit 0 ; fi
if [ ! -f $lfile ] ; then echo "locus file <$lfile> not found, stopping" ; exit 0 ; fi
if [ ! $indir ] ; then echo "input directory not specified (-d option), stopping" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "input directory <$indir> not found, stopping" ; exit 0 ; fi
if [ ! $samsuf ] ; then echo "suffix of sample subfolders not specified (-f option), setting to ''" ; samsuf="" ; fi
if [ ! $suffix ] ; then echo "suffix of locus subfolders not specified (-x option), setting to ''" ; suffix="" ; fi
 
## Additional arguments
nbcheck=$(wc -l ${sfile})
pathtocontigs="dipspades/consensus_contigs.fasta"

if [ -f nb.${suffix}.txt ] ; then /bin/rm -f nb.${suffix}.txt ; fi
if [ -f nb.assemblies.txt ] ; then /bin/rm -f nb.assemblies.txt ; fi
for i in $(cat $sfile | head -n $nbcheck)
do
	nbas=$(ls -1d ${indir}/${i}${samsuf}/*${suffix} | wc -l)
	ncon=$(ls -1 ${indir}/${i}${samsuf}/*${suffix}/${pathtocontigs} | wc -l)

	echo -e "${i}\t${nbas}" >> nb.${suffix}.txt
	echo -e "${i}\t${ncon}" >> nb.assemblies.txt

	echo "${i}: ${ncon} / ${nbas} loci assembled"
done

## Number of loci
nloci=$(wc -l < ${lfile})

## Number of assembly numbers
nnb=$(cut -f2 nb.${suffix}.txt | uniq | wc -l)
if [ "$nnb" -gt 1 ] ; then echo "found different assembly numbers!" ; cut -f2 -d ' ' nb.${suffix}.txt | uniq ; echo ; echo "number of loci: ${nloci}" ; exit 0 ; fi 

## Number of assemblies
nas=$(cut -f2 nb.${suffix}.txt | uniq)
if [ "$nas" -eq "$nloci" ] ; then echo "number of assemblies matches number of loci <${nloci}>" ; else echo "number of assemblies <${nas}> differs from number of loci <${nloci}>!" ; fi
