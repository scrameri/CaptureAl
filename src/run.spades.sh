#!/bin/bash

# best to start this from a local scratch

## Usage: run.spades.sh -s <samples.txt> -r <readpairsdir> -t <threads>

## Arguments:
# -s sample file (will look for extracted reads inside ${sample}.targets)
# -r absolute path to extracted read pairs
# -t number of threads: the parallelization is over loci, not over individuals (these are assembled one after the other)

## Needs
# submit-commands-var.pl

## Manuals
# http://cab.spbu.ru/files/release3.6.0/manual.html

## Define arguments
while getopts s:r:t: opts
do
        case "${opts}"
        in
				s) sfile=${OPTARG};;
                r) readpairsdir=${OPTARG};;
                t) threads=${OPTARG};;

    	 esac
done


## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! $readpairsdir ] ; then echo "absolute path to extracted read pairs (-r option) not specified, stopping." ; exit 0 ; fi
if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, using -t 15." ; threads=15 ; fi
if [ "$threads" -gt 30 ] ; then echo "number of threads (-t option) must be between 1 and 30, setting -t 2." ; threads=2 ; fi


## Create results directory
results=$(basename $readpairsdir)
mkdir $results
cp $sfile $results ; cp $readpairsdir/regions $results
cd $results


## Create log file
echo "=================================================================" > spades.log
echo "=========================== SPADES LOG ==========================" >> spades.log
echo "=================================================================" >> spades.log
echo " " >> spades.log
echo "Starting time:                $(zdump MEC)" >> spades.log
echo "Sample file:                  $sfile" >> spades.log
echo "Read pairs used:              ${readpairsdir}" >> spades.log
echo "Number of threads used:       ${threads}" >> spades.log


## Loop over spades samples	
for sample in $(cat ${sfile})
  do
	# create spades commands
	readlocation="${readpairsdir}/${sample}.targets"
	samplebase=$(basename ${sample})
	spadesdir="${samplebase}.spades"
	mkdir ${spadesdir}
	cd ${spadesdir}
	/bin/rm -rf ${sample}.spades.commands.txt
	
	## Loop over targets
	for fastq1 in $(ls -1 ${readlocation}/extracted_reads_*trim1*ids.fastq)
	do
	   
       	readbase1=`basename $fastq1 .fastq`
		readbase2=`echo $readbase1 |perl -p -e 's/trim1/trim2/' `
		readbase=`echo $readbase1 |perl -p -e 's/trim1.//' `
		
		## cov-cutoff auto + careful
        cmd="spades.py -1 \"${readlocation}/${readbase1}.fastq\" -2  \"${readlocation}/${readbase2}.fastq\"  --careful --cov-cutoff auto  --only-assembler  -t 1   --disable-gzip-output  -o \"${readbase}.spades\" > \"${readbase}.spades.log\" 2>\"${readbase}.spades.err\" "
        
        echo $cmd >> ${sample}.spades.commands.txt
		 	
	done
      
	## Run spades jobs:
	submit-commands-var.pl  $threads  ${sample}.spades.commands.txt
	cd ../
	
  done


## Sample Finish Time
echo "Finish time:                  $(zdump MEC)" >> spades.log
echo " " >> spades.log


## Finish
cd ../
echo 
echo "All samples processed."
echo
