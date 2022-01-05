#!/bin/bash

## Usage: get.coverage.stats.sh -d <mappingdir> -s <samples.txt> -Q <mapping quality> -t <threads>

## Needs: get.coverage.stats.R, samtools, bedtools

## Define arguments
while getopts d:s:Q:t: opts
do
        case "${opts}"
        in
        	d) in=${OPTARG};;
        	s) sfile=${OPTARG};;
        	Q) Q=${OPTARG};;
        	t) threads=${OPTARG};;
    	esac
done

## Check arguments
if [ ! $in ] ; then echo "mapping directory (-d option) not provided, using current directory" ; in=$(pwd) ; fi
if [ ! -d $in ] ; then echo "mapping directory <$in> not found, stopping" ; exit 0 ; fi
if [ ! $sfile ] ; then echo "sample file (-s option) not provided, assuming <samples.txt>" ; sfile="samples.txt" ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! $Q ] ; then echo "mapping quality threshold (-Q option) not specified (in .bam filename), stopping" ; exit 0 ; fi
if [ ! $threads ] ; then echo "number of threads (-t option) not specified, setting to 2" ; threads=2 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping." ; exit 0 ; fi
ref=$(echo ${in}/*.fasta)
if [ ! -f ${ref} ] ; then echo "reference fasta file not found in mapping directory ${in}, stopping" ; exit 0 ; fi

## Get mapping stats
if [ -f ${in}/mapping.stats.Q${Q}.txt ]
then
	echo ; echo "mapping.stats.Q${Q}.txt already exists, moving to mapping.stats.Q${Q}.txt.bak" ; echo
	mv ${in}/mapping.stats.Q${Q}.txt ${in}/mapping.stats.Q${Q}.txt.bak
fi

echo ; echo "collecting mapping stats..." ; echo
getMappingStats() {
	name=$1
	echo $name
	bamfile="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam"
	flagstats="${in}/statsQ20/${name}.flagstats.txt"
	ofile="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.mapstats.txt"

	# header	
	echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${ofile}
	
	# from samtools view: get number of unique read IDs (pairs + singletons)
	nid=$(samtools view ${bamfile} | cut -f1 | sort | uniq | wc -l)
	
	# from flagstats: get number of mapped reads (fwd + rev)
	tot=$(cat ${flagstats} | grep "in total" | cut -f1 -d' ')
	pp=$(cat ${flagstats} | grep "properly paired" | cut -f1 -d' ')
	si=$(cat ${flagstats} | grep "singletons" | cut -f1 -d' ')
	percpp=$(awk -v a="$tot" -v b="$pp" 'BEGIN {print b/a*100}')
	
	echo -e "${name}\t${nid}\t${tot}\t${pp}\t${si}\t${percpp}"  >> ${ofile}
	
}
export Q=${Q}
export -f getMappingStats
cat ${sfile} | parallel -j ${threads} getMappingStats

# collect
echo -e "sample\tread-ids\ttotal-reads\tproperly-paired\tsingletons\t%properly-paired" > ${in}/mapping.stats.Q${Q}.txt
for name in $(cat ${sfile})
do
	infile="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.mapstats.txt"
	cat $infile | tail -n+2 | head -n1 >> ${in}/mapping.stats.Q${Q}.txt
done


## Calculating coverages
echo ; echo "calculating coverages..." ; echo
doCovcalc() {
	sample=$1
	echo $name
	bamfile="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam"
     	
	## Calculate mapped length and average coverage per region
	# by default, only considers regions with properly-paired reads (argument 3)
	# by default, calculates coverages region by region to save memory in parallel mode (argument 4)	
	get.coverage.stats.R ${bamfile} ${ref}
}
export ref=${ref}
export Q=${Q}
export -f doCovcalc
cat ${sfile} | parallel -j ${threads} doCovcalc


## Get coverage stats above threshold
echo ; echo "calculating coverages above thresholds..." ; echo
doCount() {
	sample=$1
	echo ${name}
	ifile="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.coverage.txt"
	ofile="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.covtab.txt"
	/bin/rm -rf ${ofile}
	
	for i in 1 2 3 4 5 6 7 8 9 10 15 20 30 40 50 100 200 300 400 500 600 700 800 900 1000 10000
	do 
		echo -n -e  "coverage >= $i\t" >> ${ofile}
		cat ${ifile} | tail -n +2 | awk '{if($4 >='$i' ) print $0}' |wc -l >> ${ofile}
	done
}
export Q=$Q
export -f doCount
cat ${sfile} | parallel -j ${threads} doCount


## Collect coverage stats
if [ -f ${in}/coverage.stats.Q${Q}.txt ]
then
	echo "coverage.stats.Q${Q}.txt already exists, moving to coverage.stats.Q${Q}.txt.bak" ; echo
	mv ${in}/coverage.stats.Q${Q}.txt ${in}/coverage.stats.Q${Q}.txt.bak
fi

echo ; echo "collecting coverage stats..." ; echo
echo -ne "sample\t1\t10\t20\t30\t40\t50\t100\t200\t300\t400\t500\t1000\n" > ${in}/coverage.stats.Q${Q}.txt
for sample in  $(cat ${sfile})
do  
	infile="${in}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.covtab.txt"
	echo -ne "${name}\t" >> ${in}/coverage.stats.Q${Q}.txt
	for c in 1 10 20 30 40 50 100 200 300 400 500 1000 
	do 
		grep -w "coverage >= $c" ${infile} | cut -f2 | tr "\n" "\t"  >> ${in}/coverage.stats.Q${Q}.txt
	done 
	echo >> ${in}/coverage.stats.Q${Q}.txt
done


## Finish
touch ${in}/mapping.stats.Q${Q}.txt
touch ${in}/coverage.stats.Q${Q}.txt
echo ; echo "All samples processed" ; echo
