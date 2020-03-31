#!/bin/bash

## Usage: get.coverage.stats.sh -s <samples.txt> -Q <mapping quality> -t <threads>

## Needs: get.coverage.stats.R, samtools

## Define arguments
while getopts s:Q:t: opts
do
        case "${opts}"
        in
        	s) sfile=${OPTARG};;
        	Q) mapQ=${OPTARG};;
        	t) threads=${OPTARG};;
    	esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample file (-s option) not provided, stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! $mapQ ] ; then echo "mapping quality threshold (-Q option) not specified (in .bam filename), stopping" ; exit 0 ; fi
if [ ! $threads ] ; then echo "number of threads (-t option) not specified, setting to 8" ; threads=8 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping." ; exit 0 ; fi


## Collect mapping stats
if [ -f mapping.stats.txt ]
then
	echo ; echo "mapping.stats.txt already exists, moving to mapping.stats.txt.bak" ; echo
	mv mapping.stats.txt mapping.stats.txt.bak
fi

echo "\ncollecting mapping stats..." ; echo
echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > mapping.stats.txt
mkdir mapping.stats.tmp
getMappingStats() {
	sample=$1
	echo $sample
	bamfile="${sample}/${sample}.bwa-mem.mapped.Q${mapQ}.sorted.bam"
	flagstats="${sample}/${sample}.bwa-mem.mapped.Q${mapQ}.flagstats"
	
	echo -en "$sample\t"  >> mapping.stats.tmp/${sample}.tmp
	
	# from samtools view: get number of unique read IDs (pairs + singletons)
	samtools view ${bamfile} | cut -f1 | sort | uniq | wc -l | tr "\n" "\t" >> mapping.stats.tmp/${sample}.tmp
	
	# from flagstats: get number of mapped reads (fwd + rev)
	# samtools view ${bamfile} | cut -f1 | wc -l | tr "\n" "\t" >> mapping.stats.tmp/${sample}.tmp	
	cat ${flagstats} | grep -P "in total|properly paired|singletons" | cut -f1 -d "+" | tr "\n" "\t" | awk '{print $1 "\t" $2 "\t" $3 "\t" $2/$1*100}'  >> mapping.stats.tmp/${sample}.tmp
	
}
export mapQ=${mapQ}
export -f getMappingStats
cat ${sfile} | parallel -j ${threads} getMappingStats


# collect
cat mapping.stats.tmp/*tmp >> mapping.stats.txt
rm -r mapping.stats.tmp

## Calculating coverages
doCovcalc() {
	sample=$1
	echo $sample
	cd ${sample}
	ref=$(echo *.fasta)
	bamfile=${sample}.bwa-mem.mapped.Q${mapQ}.sorted.bam
     	
	## Calculate mapped length and average coverage per region
	# by default, only considers regions with properly-paired reads (argument 3)
	# by default, calculates coverages region by region to save memory in parallel mode (argument 4)	
	get.coverage.stats.R ${bamfile} ${ref}
	cd ..
}
export mapQ=${mapQ}
export -f doCovcalc
echo "\ncalculating coverages..." ; echo
cat ${sfile} | parallel -j ${threads} doCovcalc


## Get coverage stats per category
doCount() {

	sample=$1
	echo $sample
	cd $sample
	ifile="${sample}.bwa-mem.mapped.Q${mapQ}.sorted.coverage.txt"
	ofile="${sample}.coverage.above.threshold.txt"
	/bin/rm -rf ${ofile}
	
	for i in 1 2 3 4 5 6 7 8 9 10 15 20 30 40 50 100 200 300 400 500 600 700 800 900 1000 10000
	do 
		echo -n -e  "coverage >= $i\t" >> ${ofile}
		cat ${ifile} | tail -n +2 | awk '{if($4 >='$i' ) print $0}' |wc -l >> ${ofile}
	done
	cd ..
}
export -f doCount
echo "\ncalculating coverages above thresholds..." ; echo
cat ${sfile} | parallel -j ${threads} doCount


## Collect coverage stats
if [ -f coverage.stats.txt ]
then
	echo "coverage.stats.txt already exists, moving to coverage.stats.txt.bak" ; echo
	mv coverage.stats.txt coverage.stats.txt.bak
fi

echo "\ncollecting coverage stats..." ; echo
echo -ne "sample\t1\t10\t20\t30\t40\t50\t100\t200\t300\t400\t500\t1000\n" > coverage.stats.txt
for sample in  $(cat $sfile)
do  
	infile="${sample}/${sample}.coverage.above.threshold.txt"
	echo -ne "${sample}\t" >> coverage.stats.txt
	for c in 1 10 20 30 40 50 100 200 300 400 500 1000 
	do 
		grep -w "coverage >= $c" ${infile} | cut -f2 | tr "\n" "\t"  >> coverage.stats.txt
	done 
	echo >> coverage.stats.txt
done


## Finish
touch mapping.stats.txt
touch coverage.stats.txt
echo ; echo "All samples processed" ; echo
