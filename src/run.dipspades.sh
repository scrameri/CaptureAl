#!/bin/bash

# best to start this from a local scratch

## Usage: run.dipspades.sh -s <samples.txt> -r <readpairsdir> -t <threads>

## Arguments:
# -s sample file (will look for extracted reads inside ${name}.targets)
# -r absolute path to extracted read pairs
# -t number of threads: the parallelization is over reference sequences, not over individuals (these are assembled one after the other)

## Needs
# submit-commands-var.pl

## Define arguments
while getopts s:r:t: opts
do
        case "${opts}"
        in
				s) sfile=${OPTARG};;
                r) in=${OPTARG};;
                t) threads=${OPTARG};;

    	 esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! $in ] ; then echo "absolute path to extracted read pairs (-r option) not specified, stopping." ; exit 0 ; fi
if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, using -t 15." ; threads=15 ; fi
if [ "$threads" -gt 30 ] ; then echo "number of threads (-t option) must be between 1 and 30, setting -t 15." ; threads=15 ; fi

## Create output directory
out=$(basename ${in})
mkdir ${out}
cp ${sfile} ${out} ; cp ${in}/loci.txt ${out}
cd ${out}

## Create log file
logfile="dipspades.log"

echo "=================================================================" > ${logfile}
echo "========================= DIPSPADES LOG =========================" >> ${logfile}
echo "=================================================================" >> ${logfile}
echo " " >> ${logfile}
echo "Starting time:                $(zdump CET)" >> ${logfile}
echo "Sample file:                  $sfile" >> ${logfile}
echo "Read pairs used:              ${in}" >> ${logfile}
echo "Number of threads used:       ${threads}" >> ${logfile}

## Loop over dipspades samples	
for name in $(cat ${sfile})
  do
	# create dipspades commands
	reads="${in}/${name}.targets"
	odir="${name}.dipspades"
	mkdir ${odir}
	cd ${odir}
	/bin/rm -rf ${name}.dipspades.commands.txt
	
	## Loop over targets
	for fastq1 in $(ls -1 ${reads}/extracted_reads_*trim1*ids.fastq)
	do
	   
       	readbase1=`basename $fastq1 .fastq`
		readbase2=`echo $readbase1 |perl -p -e 's/trim1/trim2/' `
		readbase=`echo $readbase1 |perl -p -e 's/trim1.//' `
		
		## cov-cutoff auto + careful
        cmd="~/bin/SPAdes-3.6.0-Linux/bin/dipspades.py -1 \"${reads}/${readbase1}.fastq\" -2  \"${reads}/${readbase2}.fastq\"  --careful --cov-cutoff auto  --only-assembler  -t 1   --disable-gzip-output  -o \"${readbase}.spades\" > \"${readbase}.dipspades.log\" 2>\"${readbase}.dipspades.err\" ; rm -rf \"${readbase}.spades/spades\" "
        
        ## cov-cutoff off + no careful
        #cmd="~/bin/SPAdes-3.6.0-Linux/bin/dipspades.py -1 \"${reads}/${readbase1}.fastq\" -2  \"${reads}/${readbase2}.fastq\" --cov-cutoff off  --only-assembler  -t 1   --disable-gzip-output  -o \"${readbase}.spades\" > \"${readbase}.dipspades.log\" 2>\"${readbase}.dipspades.err\" ; rm -rf \"${readbase}.spades/spades\" "
        echo $cmd >> ${name}.dipspades.commands.txt
		 	
	done
      
	## Run dipspades jobs:
	submit-commands-var.pl  $threads  ${name}.dipspades.commands.txt
	cd ../
	
  done

## Sample Finish Time
echo "Finish time:                  $(zdump CET)" >> ${logfile}
echo " " >> ${logfile}

## Finish
cd ../
echo 
echo "All samples processed."
echo
