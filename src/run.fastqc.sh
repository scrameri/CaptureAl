#!/bin/bash

## Define arguments
while getopts s:a:d:o:t: opts
do
        case "${opts}"
        in
                s) sfile=${OPTARG};;
                a) adapters=${OPTARG};;
                d) in=${OPTARG};;
                o) out=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done


# check input
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping." ; exit 0 ; fi

if [ ! $adapters ] ; then echo "adapter file (-a option) not specified, stopping" ; exit 0 ; fi
if [ ! -f $adapters ] ; then echo "adapter file <$adapters> not found, stopping." ; exit 0 ; fi

if [ ! $in  ] ; then echo "path to reads (-d option) not specified, using current directory" ; in=$currentdir ; fi
if [ ! $out  ] ; then echo "path to output directory (-o option) not specified, using current directory" ; out=$currentdir ; fi

if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, using -t 16." ; threads=16 ; fi
if [ "$threads" -gt 30 ] ; then echo "number of threads (-t option) too large, setting -t 16." ; threads=16 ; fi


# create output directory if needed
if [ ! -d $(dirname ${out}) ] ; then mkdir $(dirname ${out}) ; fi
if [ ! -d ${out} ] ; then mkdir ${out} ; fi


# define function
doFASTQC() {
	
	name=$1

	reads1=$(ls -1d ${in}/${name}[_.]R1[_.]*gz)
	reads2=$(ls -1d ${in}/${name}[_.]R2[_.]*gz)

	trim1=$(ls -1d ${in}/${name}.trim1.fastq.gz)
	trim2=$(ls -1d ${in}/${name}.trim2.fastq.gz)

	unp1=$(ls -1d ${in}/${name}.U1.fastq.gz)
	unp2=$(ls -1d ${in}/${name}.U2.fastq.gz)

	# run fastqc
	for file in $reads1 $reads2 $trim1 $trim2 $unp1 $unp2
	do
	
		if [ -f $file ]
		then
			fastqc $file -t 1 -a ${adapters} -o ${out}
		fi
	
	done
	
}

# export
export in=$in
export out=$out
export adapters=$adapters

cat ${sfile} | parallel -j $threads doFASTQC
