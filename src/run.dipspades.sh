#!/bin/bash

# best to start this from a local scratch

## Usage: run.dipspades.sh -s <samples.txt> -r <readpairsdir> -t <threads>

## Needs
# submit-commands-var.pl

## Arguments:
# -s sample file (will look for extracted reads inside ${sample}.targets)
# -r absolute path to extracted read pairs
# -t number of threads: the parallelization is over reference sequences, not over individuals (these are assembled one after the other)

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
if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, setting to -t 4." ; threads=4 ; fi


## Create results directory
results=$(basename $readpairsdir)
mkdir $results
cp $sfile $results
cd $results


## Create log file
echo "=================================================================" > dipspades.log
echo "========================= DIPSPADES LOG =========================" >> dipspades.log
echo "=================================================================" >> dipspades.log
echo " " >> dipspades.log
echo "Starting time:                $(zdump MEC)" >> dipspades.log
echo "Sample file:                  $sfile" >> dipspades.log
echo "Read pairs used:              ${readpairsdir}" >> dipspades.log
echo "Number of threads used:       ${threads}" >> dipspades.log


## Loop over dipspades samples	
for sample in $(cat ${sfile})
  do
	# create dipspades commands
	readlocation="${readpairsdir}/${sample}.targets"
	samplebase=$(basename ${sample})
	dipspadesdir="${samplebase}.dipspades"
	mkdir ${dipspadesdir}
	cd ${dipspadesdir}
	/bin/rm -rf ${sample}.dipspades.commands.txt
	
	## Loop over targets
	for fastq1 in $(ls -1 ${readlocation}/extracted_reads_*trim1*ids.fastq)
	do
	   
       	readbase1=`basename $fastq1 .fastq`
	readbase2=`echo $readbase1 |perl -p -e 's/trim1/trim2/' `
	readbase=`echo $readbase1 |perl -p -e 's/trim1.//' `
		
	## cov-cutoff auto + careful
        cmd="~/bin/SPAdes-3.6.0-Linux/bin/dipspades.py -1 \"${readlocation}/${readbase1}.fastq\" -2  \"${readlocation}/${readbase2}.fastq\"  --careful --cov-cutoff auto  --only-assembler  -t 1   --disable-gzip-output  -o \"${readbase}.spades\" > \"${readbase}.dipspades.log\" 2>\"${readbase}.dipspades.err\" ; rm -rf \"${readbase}.spades/spades\" "
        
        ## cov-cutoff off + no careful
        #cmd="~/bin/SPAdes-3.6.0-Linux/bin/dipspades.py -1 \"${readlocation}/${readbase1}.fastq\" -2  \"${readlocation}/${readbase2}.fastq\" --cov-cutoff off  --only-assembler  -t 1   --disable-gzip-output  -o \"${readbase}.spades\" > \"${readbase}.dipspades.log\" 2>\"${readbase}.dipspades.err\" ; rm -rf \"${readbase}.spades/spades\" "
        echo $cmd >> ${sample}.dipspades.commands.txt
		 	
	done
      
	## Run dipspades jobs:
	submit-commands-var.pl  $threads  ${sample}.dipspades.commands.txt
	cd ../
	
  done


## Sample Finish Time
echo "Finish time:                  $(zdump MEC)" >> dipspades.log
echo " " >> dipspades.log


## Finish
cd ../
echo 
echo "All samples processed."
echo


