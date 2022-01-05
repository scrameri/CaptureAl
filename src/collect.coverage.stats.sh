#!/bin/bash

## Usage: collect.coverage.stats.sh -s <samples.txt> -Q <mapping quality> -m <maximum coverage> -d <mappingdir>

## Needs: R
module --silent load r 1>/dev/null 2>/dev/null

## Define arguments
while getopts d:s:Q:m: opts
do
        case "${opts}"
        in
        	d) in=${OPTARG};;
        	s) sfile=${OPTARG};;
        	Q) Q=${OPTARG};;
        	m) maxcov=${OPTARG};;
    	esac
done

## Arguments
if [ ! ${sfile} ] ; then echo "sample file (-s option) not provided, assuming <samples.txt>" ; sfile="samples.txt" ; fi
if [ ! -f ${sfile} ] ; then echo "sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! ${Q} ] ; then echo "mapping quality threshold (-Q option) not specified (see .bam filename), assuming Q=<10>" ; Q=10 ; fi
if [ ! ${maxcov} ] ; then echo "maximum coverage threshold (-m option) not specified (see .bam filename), assuming maxcov=<500>" ; maxcov=500 ; fi
if [ ! ${in} ] ; then echo "mapping directory (-d option) not provided, using current directory" ; in=$(pwd) ; fi
if [ ! -d ${in} ] ; then echo "mapping directory <${in}> not found, stopping" ; exit 0 ; fi

## Collect mapping stats
echo ; echo "collecting mapping stats..."

echo -e "sample\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${in}/mapping.stats.txt
echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${in}/mapping.stats.Q${Q}.txt
echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${in}/mapping.stats.Q${Q}.nodup.cov${maxcov}.txt

for name in $(cat ${sfile})
do
	mstatfile="${in}/${name}/${name}.mapstats.txt"
	mstatfileQ="${in}/${name}/${name}.Q${Q}.mapstats.txt"
	mstatfileNODUP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.mapstats.txt"

	if [ -f ${mstatfile} ]
	then
		cat ${mstatfile} | tail -n+2 | head -n1 >> ${in}/mapping.stats.txt
	else
		echo "<${mstatfile}> not found for ${name}"
	fi

	if [ -f ${mstatfileQ} ]
	then
		cat ${mstatfileQ} | tail -n+2 | head -n1 >> ${in}/mapping.stats.Q${Q}.txt
	else
		echo "<${mstatfileQ}> not found for ${name}"
	fi
	
	if [ -f ${mstatfileNODUP} ]
	then
		cat ${mstatfileNODUP} | tail -n+2 | head -n1 >> ${in}/mapping.stats.Q${Q}.nodup.cov${maxcov}.txt
	else
		echo "<${mstatfileNODUP}> not found for ${name}"
	fi
done


## Collect coverage stats
echo ; echo "collecting coverage above threshold stats..."

echo -ne "sample\t1\t10\t20\t30\t40\t50\t100\t200\t300\t400\t500\t1000\n" > ${in}/coverage.stats.Q${Q}.txt
echo -ne "sample\t1\t10\t20\t30\t40\t50\t100\t200\t300\t400\t500\t1000\n" > ${in}/coverage.stats.Q${Q}.nodup.cov${maxcov}.txt

for name in  $(cat ${sfile})
do  
	covtabQ="${in}/${name}/${name}.Q${Q}.covtab.txt"
	covtabNODUP="${in}/${name}/${name}.Q${Q}.nodup.cov${maxcov}.covtab.txt"
	
	if [ -f ${covtabQ} ]
	then
		echo -ne "${name}\t" >> ${in}/coverage.stats.Q${Q}.txt
		for c in 1 10 20 30 40 50 100 200 300 400 500 1000 
		do 
			grep -w "coverage >= $c" ${covtabQ} | cut -f2 | tr "\n" "\t"  >> ${in}/coverage.stats.Q${Q}.txt
		done
		echo >> ${in}/coverage.stats.Q${Q}.txt
	else
		echo "<${covtabQ}> not found for ${name}"
	fi

	
	if [ -f ${covtabNODUP} ]
	then
		echo -ne "${name}\t" >> ${in}/coverage.stats.Q${Q}.nodup.cov${maxcov}.txt
		for c in 1 10 20 30 40 50 100 200 300 400 500 1000 
		do 
			grep -w "coverage >= $c" ${covtabNODUP} | cut -f2 | tr "\n" "\t"  >> ${in}/coverage.stats.Q${Q}.nodup.cov${maxcov}.txt
		done
		echo >> ${in}/coverage.stats.Q${Q}.nodup.cov${maxcov}.txt
	else
		echo "<${covtabNODUP}> not found for ${name}"
	fi
done

## Collect coverage stats
echo
collect.coverage.stats.R ${sfile} ${Q} ${maxcov} ${in}
