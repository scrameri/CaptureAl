#!/bin/bash

################################################################################
### RUN 'bwa mem' on multiple samples automatically and in multi-thread mode ###

## USAGE: 
# run.bwamem.sh -s <samples.txt> -r <file> -e <string(s) separated by ','> -T <positive integer> -Q <positive integer> -d <path to input reads> -o <output directory> -t <positive integer>

## DETAILS:
#   sample file (-s option) must contain the basenames of all samples to be mapped (basenames = filenames without file extensions, e.g. 'SH598_S16')
#   sample file extensions (-e option) of files 'SH598_S16.trim1.fastq.gz' and 'SH598_S16.trim2.fastq.gz' might be '.trim1.fastq.gz,.trim2.fastq.gz'
#   -> only list one basename per sample for paired data (where the two files have the same basename, but different file extensions, see -e option)

# Required arguments:
# -r reference fasta file
# -e sample file(s) extension(s). E.g. '.fasta' or '.trim.fq.gz' [unpaired] or '.trim1.fastq,.trim2.fastq' [paired]. The program then interprets if data is unpaired or paired. Separate file extensions of file pairs with a ','.

# Optional arguments:
# -T [10]  bwa mem -T alignment score [mapped reads will go to $TMPDIR]
# -Q [20]  minimum mapping quality (fifth field / MAPQ in .sorted.bam file) [high-quality reads will go to $TMPDIR]. Only deduplicated reads will go to ${out}/${name}/*bwa-mem.sorted.Q${Q}.nodup.bam. Must be Q >= T. If above 0, will filter reads with multiple mappings.
# -d [pwd] path to input reads
# -o [pwd] path to output directory (will be created if it does not exist)
# -t [3]   number of samples processed in parallel. Can be between 1 (uses ${cpu} CPU cores in total) and 6 (uses 6*${cpu} CPU cores in total, cpu=4 by default)


## Needs: bwa, sambamba, java, python, picard-tools, r

## OUTPUT: folders <sample> containing .log, .sorted.Q${Q}.nodup.bam(.bai) and stat folders <stats> (flagstats of mapped reads), <statsQ${Q} (flagstats of high-quality reads), <stats_dup> (PCR duplication stats)

## Authors
# Simon Crameri, ETHZ
# Stefan Zoller, GDC

################################################################################

## Define arguments
while getopts s:r:T:Q:e:d:o:t: opts
do
        case "${opts}"
        in
                s) sfile=${OPTARG};;
                r) fas=${OPTARG};;
                T) T=${OPTARG};;
                Q) Q=${OPTARG};;
                e) ext=${OPTARG};;
                d) in=${OPTARG};;
                o) out=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done

# load modules
module load gdc gcc/4.8.2 bwa/0.7.17 sambamba/0.8.0 java/1.8.0_101 python/3.6.1 java/1.8.0_73 picard-tools/2.23.8 r/3.1.2

# check required files and arguments and set further arguments 
currentdir=$(pwd)
ext1=$(echo ${ext} | cut -d "," -f1)
ext2=$(echo ${ext} | cut -d "," -f2)

fasext=${fas##*.} # $fas file extension
fasline=$(head -n1 ${fas}) # first line in $fas 
faschr=$(echo "${fasline:0:1}") # should be a '>'
cpu=4 # parallel processing of bwa mem

# check input
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping." ; exit 0 ; fi
if [ ! $fas ] ; then echo "reference fasta file not specified (-r option), stopping." ; exit 0 ; fi
if [[ $fasext = "fasta" ]] ; then echo "" ; else echo " reference fasta file (-r option) does not have .fasta extension, stopping." ; exit 0 ; fi
if [[ $faschr = ">" ]] ; then echo "" ; else echo "reference fasta (-r option) seems not to be in .fasta format, stopping." ; exit 0 ; fi
if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, using -t 3." ; threads=3 ; fi
if [ "$threads" -gt 6 ] ; then echo "number of threads (-t option) must be between 1 (uses $cpu CPU cores in total) and 6 (uses 6*$cpu CPU cores in total), setting -t 3." ; threads=3 ; fi
if [ ! $T ] ; then echo "BWA MEM quality parameter (-T option) not specified, using -T 10." ; T=10 ; fi
if [ ! $Q ] ; then echo "mapping quality parameter (-Q option) not specified, using -Q 20." ; Q=20 ; fi
if [ ! $ext ] ; then echo "sample file extensions (-e option) not specified, stopping." ; exit 0 ; fi
if [ ! $in  ] ; then echo "path to reads (-d option) not specified, using current directory" ; in=$currentdir ; fi
if [ ! $out  ] ; then echo "path to output directory (-o option) not specified, using current directory" ; out=$currentdir ; fi
if [ ! -d $in ] ; then echo "path to reads <$in> not found, stopping" ; exit 0 ; fi
if [ $ext1 = $ext2 ] ; then mode='unpaired' ; else mode='paired' ; fi
if [[ $(tail -n 1 ${sfile}) = "" ]] ; then echo "found empty last line in ${sfile}, please check your ${sfile} file, stopping." ; exit 0 ; fi

# create output directories if needed
if [ ! -d $(dirname $out) ] ; then mkdir $(dirname $out) ; fi
if [ ! -d ${out} ] ; then mkdir ${out} ; fi
if [ ! -d ${out}/stats ] ; then mkdir ${out}/stats ; fi
if [ ! -d ${out}/statsQ${Q} ] ; then mkdir ${out}/statsQ${Q} ; fi
if [ ! -d ${out}/stats_dup ] ; then mkdir ${out}/stats_dup ; fi
if [ ! -f ${out}/$(basename ${sfile}) ] ; then cp ${sfile} ${out} ; fi

# check arguments
if [ "$T" -gt "$Q" ] ; then echo "high-quality parameter (-Q) must be equal or greater than mapping quality parameter (-T), stopping." ; exit 0 ; fi

# copy reference.fasta and fix headers if needed (replace ' ' and '|' with '_')
ref="${out}/$(basename ${fas})"
if [ ! -f ${ref} ] ; then cp ${fas} ${out} ; fi
fix.fasta.headers.R $ref ' |\\|' '_' FALSE

# index reference fasta
if [ ! -f ${ref}.bwt ] ; then bwa index ${ref} ; fi

# define mapping function
doMapping() {
	name=$1
	echo "mapping ${name}"
	
	# create sample directory in output directory $out
	if [ ! -d ${out}/${name} ] ; then mkdir ${out}/${name} ; fi
	
	# create log file
	logfile=${out}/${name}/bwa.log

	echo "=================================================================" > $logfile
	echo "========================== BWA MEM LOG ==========================" >> $logfile
	echo "=================================================================" >> $logfile
	echo " " >> $logfile
	echo "Starting time:                $(zdump CET)" >> $logfile
	echo "Working directory:            $(pwd)" >> $logfile
	echo "Input directory:              ${in}" >> $logfile
	echo "Output directory:             ${out}" >> $logfile
	echo "Sample file:                  ${sfile}" >> $logfile
	echo "Reference used for BWA MEM:   ${fas}" >> $logfile
	echo "File extension:               ${ext}" >> $logfile
	echo "Mode:                         ${mode}" >> $logfile
	echo "Path to sample:               ${in}/${name}${ext1}" >> $logfile
	if [ $mode = 'paired' ] ; then echo "Path to sample [pair]:        ${in}/${name}${ext2}" >> $logfile ; fi
	echo "-T bwa mem quality parameter: $T" >> $logfile
	echo "-Q mapping quality parameter: $Q" >> $logfile

	## MAPPING: BWA MEM
	# some bwa mem options
	# t: threads
	# M: flag shorter split hits as secondary (important for picard tools)
	# w: band width, essentially gaps longer than this are not found (default 100)
	# T: dont output with score lower than T (default 30)
	# B: mismatch penalty (default 4)
	# E: gap extension penalty (default 1)
	
	# map reads using bwa mem
	if [ $mode = 'paired' ]
	then

		bwa mem -t ${cpu} -M -T ${T} -R "@RG\tID:$name\tSM:$name" $ref ${in}/${name}${ext1} ${in}/${name}${ext2} > ${TMPDIR}/${name}.bwa-mem.sam 2> ${out}/${name}/${name}.bwa-mem.log

	elif [ $mode = 'unpaired' ]
	then
	
		bwa mem -t ${cpu} -M -T ${T} -R "@RG\tID:$name\tSM:$name" $ref ${in}/${name}${ext1} > ${TMPDIR}/${name}.bwa-mem.sam 2> ${out}/${name}/${name}.bwa-mem.log
	
	fi

	# convert .sam to sorted.bam
	sambamba view -t ${cpu} -S ${TMPDIR}/${name}.bwa-mem.sam -f bam -o /dev/stdout |sambamba sort /dev/stdin -o /dev/stdout -t ${cpu} -l 0 -m 6GB --tmpdir ${TMPDIR} > ${TMPDIR}/${name}.bwa-mem.sorted.bam

	# mapping flagstats
	sambamba flagstat ${TMPDIR}/${name}.bwa-mem.sorted.bam > ${out}/stats/${name}.flagstats.txt

	# remove reads with low mapping quality
	sambamba view -F "mapping_quality >= ${Q}" ${TMPDIR}/${name}.bwa-mem.sorted.bam -o ${TMPDIR}/${name}.bwa-mem.sorted.Q${Q}.bam -t ${cpu} -f bam 

	# remove PCR duplicates and write to sample directory
	#picard MarkDuplicates TMP_DIR=${TMPDIR} INPUT=${TMPDIR}/${name}.bwa-mem.sorted.Q${Q}.bam OUTPUT=${out}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam -M ${out}/stats_dup/${name}.dupstats.txt -VALIDATION_STRINGENCY LENIENT -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 512 -REMOVE_DUPLICATES true
	picard MarkDuplicates -I ${TMPDIR}/${name}.bwa-mem.sorted.Q${Q}.bam -O ${out}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam -M ${out}/stats_dup/${name}.dupstats.txt -TMP_DIR ${TMPDIR} -VALIDATION_STRINGENCY LENIENT -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 512 -REMOVE_DUPLICATES true

	# create index file
	sambamba index ${out}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam

	# final mapping flagstats
	sambamba flagstat ${out}/${name}/${name}.bwa-mem.sorted.Q${Q}.nodup.bam > ${out}/statsQ${Q}/${name}.flagstats.txt
        
	# log finish Time
	echo "Finish time:                  $(zdump CET)" >> $logfile
	echo " " >> $logfile

	touch ${out}/${name}

	
}
    
## Export variables (they are not available inside the doMapping() function otherwise)
export in=$in
export out=$out
export sfile=$sfile
export fas=$fas
export ext1=$ext1
export ext2=$ext2
export ref=$ref
export ext1=$ext1
export ext2=$ext2
export T=$T
export Q=$Q
export cpu=$cpu
export ext=$ext
export mode=$mode
export -f doMapping

## DO MAPPING
echo
echo "Start mapping..."
echo
cat ${sfile} | parallel -j $threads doMapping

## Finish
echo 
echo "All samples processed."
echo
