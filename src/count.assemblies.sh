#!/bin/bash

## Usage: count.assemblies.sh -d <assemblydir> -s <samples.txt> -o <output file> -c <region file>

## Details
# -d	path to directory with sample subdirectories with .tar.gz assemblies (default: pwd)
# -s	path to file with sample basenames (default: samples.txt)
# -o	name of output file (required)
# -c	path to file with extracted_reads_SAMPLE.${region}.spades.tar.gz names (default: all .spades.tar.gz in sample subdirectory)

## Needs: R
module load r

## Define arguments
while getopts d:s:o:c: opts
do
        case "${opts}"
        in
        	d) in=${OPTARG};;
        	s) sfile=${OPTARG};;
        	o) ofile=${OPTARG};;
        	c) cfile=${OPTARG};;
    	esac
done

## Check arguments
if [ ! ${in} ] ; then echo "assembly directory (-d option) not provided, using current directory" ; in=$(pwd) ; fi
if [ ! -d ${in} ] ; then echo "assembly directory <${in}> not found, stopping" ; exit 0 ; fi
if [ ! ${sfile} ] ; then echo "sample file (-s option) not provided, assuming <samples.txt>" ; sfile="samples.txt" ; fi
if [ ! -f ${sfile} ] ; then echo "sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! ${ofile} ] ; then echo "output file (-o option) not provided, stopping" ; exit 0 ; fi
if [ -f ${ofile} ] ; then echo "output file <$ofile> exists, moving to <$ofile.bak>" ; mv ${ofile} ${ofile}.bak ; fi
if [ ! ${cfile} ] ; then echo "region file / regex (-c option) not provided, assuming all regions in <${in}>" ; fi

## Additional arguments
astring="contigs.fasta" # search for this file in each assembly directory

## Count assemblies
echo
echo -e "sample\tassembled\ttotal"
echo -e "sample\tassembled\ttotal" > ${ofile}
for i in $(cat ${sfile})
do	

	if [ ! ${cfile} ]
	then
		ls -1d ${in}/${i}/*spades.tar.gz | Rscript -e 'write.table(basename(as.character(read.table("stdin")[,1])), file = "tmp.txt", quote = F, row.names = F, col.names = F)'
		cat tmp.txt |sed -e "s/${i}/SAMPLE/g" > compare.${i}.txt
		/bin/rm -f tmp.txt
		f="compare.${i}.txt"
	elif [ ! -f ${cfile} ]
	then
		f=$(echo ${cfile} | sed -e "s/SAMPLE/${i}/")
		if [ ! -f ${f} ] ; then echo "file <${f}> not found, stopping" ; exit 0 ; fi
		if [ ! -f compare.${i}.txt ] ; then echo "file <compare.${i}.txt> not found, stopping" ; exit 0 ; fi
	elif [ -f ${cfile} ]
	then
		f=${cfile}
	fi
	
	tot=$(cat ${f} | wc -l)	
	ass=$(
		for j in $(cat ${f} )
		do
			k=$(echo ${j} | sed -e "s/SAMPLE/${i}/g")
			tar -tvf ${in}/${i}/${k}
		done | grep -c ${astring}
		)
	
	echo -e "${i}\t${ass}\t${tot}"
	echo -e "${i}\t${ass}\t${tot}" >> ${ofile}

done
