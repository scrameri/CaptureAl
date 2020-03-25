#!/bin/bash

################################################################################
### RUN 'bwa mem' on multiple samples automatically and in multi-thread mode ###

## USAGE: 
# run.bwamem.sh -s <samples.txt> -r <file> -e <string(s) separated by ','> -T <positive integer> -Q <positive integer> -d <path> -t <positive integer>

## DETAILS:
#   sample file (-s option) must contain the basenames of all samples to be mapped (basenames = filenames without file extensions, might be 'SH598_S16_L001')
#   sample file extensions (-e option) of files 'SH598_S16_L001.trim1.fastq.gz' and 'SH598_S16_L001.trim2.fastq.gz' might be '.trim1.fastq.gz,.trim2.fastq.gz'
#   -> just one name per sample for paired data (where the two files have the same basename, but different file extensions, see -e option)

# Required arguments:
# -r reference fasta file
# -e sample file(s) extension(s). E.g. '.fasta' or '.trim.fq.gz' [unpaired] or '.trim1.fastq,.trim2.fastq' [paired]. The program then interprets if data is unpaired or paired. Separate file extensions of file pairs with a ','.

# Optional arguments:
# -T [10]  bwa mem -T mapping quality parameter [quality of .mapped.sorted.bam output]
# -Q [20]  high-quality portion mapping parameter, must be >= -T. Defines the portion of reads considered to have mapped with high quality [quality of .mapped.Q${Q}.sorted.bam]
# -d [pwd] /path/to/samples
# -t [3]   number of samples processed in parallel. Can be between 1 (uses 5 CPU cores in total) and 4 (uses 20 CPU cores in total)


## OUTPUT: folders <sample> containing .bam, .mapped.sorted.bam(.bai), .mapped.Q${Q}.sorted.bam(.bai), .mapped.sorted.Q${Q}.bam.html, .flagstats and .log files 

## Authors
# Simon Crameri, ETHZ
# Stefan Zoller, GDC

################################################################################

## Define arguments
while getopts s:r:T:Q:e:d:t: opts
do
        case "${opts}"
        in
                s) sfile=${OPTARG};;
                r) ref=${OPTARG};;
                T) qual=${OPTARG};;
                Q) Q=${OPTARG};;
                e) ext=${OPTARG};;
                d) samplepath=${OPTARG};;
                t) threads=${OPTARG};;
    	 esac
done


## Check required files and arguments and set further arguments 
currentdir=$(pwd)
refname=$(basename $ref .fasta)
ext1=$(echo $ext | cut -d "," -f1)
ext2=$(echo $ext | cut -d "," -f2)

rname=$(basename $ref) # $ref w/o full path
refext=${rname##*.} # $ref file extension
refline=$(head -n1 $ref) # first line in $ref 
refchr=$(echo "${refline:0:1}") # should be a '>'


echo
if [ ! $sfile ] ; then echo "sample file (-s option) not specified, stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping." ; exit 0 ; fi
if [ ! $ref ] ; then echo "reference fasta file not specified (-r option), stopping." ; exit 0 ; fi
if [[ $refext = "fasta" ]] ; then echo "" ; else echo " reference fasta file (-r option) does not have .fasta extension, stopping." ; exit 0 ; fi
if [[ $refchr = ">" ]] ; then echo "" ; else echo "reference fasta (-r option) seems not to be in .fasta format, stopping." ; exit 0 ; fi
if [ ! $threads  ] ; then echo "number of threads (-t option) not specified, using -t 3." ; threads=3 ; fi
if [ "$threads" -gt 6 ] ; then echo "number of threads (-t option) must be between 1 (uses 5 CPU cores in total) and 6 (uses 30 CPU cores in total), setting -t 3." ; threads=3 ; fi
if [ ! $qual ] ; then echo "BWA MEM quality parameter (-T option) not specified, using -T 10." ; qual=10 ; fi
if [ ! $Q ] ; then echo "mapping quality parameter (-Q option) not specified, using -Q 20." ; Q=20 ; fi
if [ ! $ext ] ; then echo "sample file extensions (-e option) not specified, stopping." ; exit 0 ; fi
if [ ! $samplepath  ] ; then echo "path to reads (-d option) not specified, setting to current directory" ; samplepath=$currentdir ; fi
if [ ! -d $samplepath ] ; then echo "path to reads <$samplepath> not found, stopping" ; exit 0 ; fi
if [ $ext1 = $ext2 ] ; then mode='unpaired' ; else mode='paired' ; fi
if [[ $(tail -n 1 ${sfile}) = "" ]] ; then echo "found empty last line in ${sfile}, please check your ${sfile} file, stopping." ; exit 0 ; fi


## Check input
# create temporary files
if [ -e tmp.headers ] ; then /bin/rm -f tmp.headers ; touch tmp.headers ; else touch tmp.headers ; fi
if [ -e tmp.badchar ] ; then /bin/rm -f tmp.badchar ; touch tmp.badchar ; else touch tmp.badchar ; fi
if [ -e tmp.folderlist ] ; then /bin/rm -f tmp.folderlist ; touch tmp.folderlist ; else touch tmp.folderlist ; fi
if [ -e tmp.check ] ; then /bin/rm -f tmp.check ; touch tmp.check ; else touch tmp.check ; fi


# check that -Q >= -T
if [ "$qual" -gt "$Q" ] ; then echo "high-quality parameter (-Q option) must be equal or greater than mapping quality parameter (-T option), stopping." ; exit 0 ; fi

# check reference fasta headers
grep '^>' $ref > tmp.headers

while read i
do 
	t1=$(echo $i | cut -d " " -f1)
	t2=$(echo $i | cut -d "|" -f1)
	if [[ ( $t1 != $i ) || ( $t2 != $i ) ]] ; then echo "stop" >> tmp.badchar ; fi
done < tmp.headers

ncheck=$(wc -l tmp.badchar | cut -d " " -f1)

if [ $ncheck != 0 ]
then	
	
	echo
	read -r -p "Found ' ' and/or '|' characters in reference fasta file headers, these may cause troubles in downstream analyses. Continue anyway? [y/n] " response
	case $response in
		[nN][oO]|[nN])
		
		/bin/rm -f tmp.headers
		/bin/rm -f tmp.badchar
		echo "Stopping." 
		exit 0
		
	esac
	
	echo
	read -r -p "Would you like to fix the reference fasta file headers (replace ' ' and/or '|' characters with '_', original file is kept)? [y/n] " response
	case $response in
		[yY][eE][sS]|[yY])
		
		for i in `cat tmp.headers`
  		do
			# Change fasta headers (old file is saved as .badnames)
			cp $ref ${refname}.badheaders.fasta
			sed -e 's/ /_/g' $ref > tmp.ref.corr
			sed -e 's/|/_/g' tmp.ref.corr > tmp.ref.corr2
			mv tmp.ref.corr2 ${ref}
			/bin/rm -f tmp.ref.corr
 		done

	esac
fi


# check existing dirs
for d in `cat ${sfile}`
do 
	sample=`basename $d`
	#echo "${sample}.on.${refname}.bwa" >> tmp.folderlist
	#if [ -e ${sample}.on.${refname}.bwa ] ; then echo "stop" >> tmp.check ; fi
	echo "${sample}" >> tmp.folderlist
	if [ -e ${sample} ] ; then echo "stop" >> tmp.check ; fi
done
ncheck=$(wc -l tmp.check | cut -d " " -f1)

if [ $ncheck != 0 ]
then
	echo
	read -r -p "One or more output folders already exist. Replace with new output folders? [y/n] " response
	case $response in
		[nN][oO]|[nN])
	
		/bin/rm -f tmp.headers
		/bin/rm -f tmp.badchar
		/bin/rm -f tmp.folderlist
		/bin/rm -f tmp.check
		echo "Stopping." 
		exit 0
	esac
fi


# remove dirs with same name as output dirs
while IFS= read -r f
do
    if [ -e $currentdir/$f ] ; then /bin/rm -rf $f ; fi
done < tmp.folderlist


# remove temporary files
/bin/rm -f tmp.check
/bin/rm -f tmp.headers
/bin/rm -f tmp.folderlist
/bin/rm -f tmp.badchar


## Build mapping index
echo
echo "Creating bwa index files"
refbase=`basename $ref` 
if [ -e $currentdir/$refbase ] ; then bwa index $refbase ; else cp $ref $currentdir ; bwa index $refbase ; fi
echo


## Make dirs and link input files
for d in $(cat ${sfile})
do
	sample=`basename $d`
	echo "making dir for $sample"
	mkdir ${sample} #.on.${refname}.bwa
	cd ${sample}  #.on.${refname}.bwa
	
	if [ -e ${samplepath}/${sample}${ext1} ]
	then 
		:
	else 
		echo "${samplepath}/${sample}${ext1} does not exist."
		echo "Check your file extensions (-e argument) or path/to/samples (-d argument)."
		exit
	fi
	
	if [ $mode = 'paired' ]
	then 
		if [ -e ${samplepath}/${sample}${ext2} ]
		then 
			:
		else 
			echo "${samplepath}/${sample}${ext2} does not exist." 
			echo "Check your file extensions (-e argument) or path/to/samples (-d argument)."
			exit
		fi
	fi
	
	ln -s ${samplepath}/${sample}${ext1} .                          	    
	if [ $mode = 'paired' ] ; then ln -s ${samplepath}/${sample}${ext2} . ; fi 
	cd ..
done


## Define mapping function
doMapping() {
	sample=$1
	echo "mapping $sample "
   	cd ${sample} #.on.${refname}.bwa
	
	ln -s ../${ref}* .
	
   	## Create LOG file
	echo "=================================================================" > bwa.log
	echo "========================== BWA MEM LOG ==========================" >> bwa.log
	echo "=================================================================" >> bwa.log
	echo " " >> bwa.log
	echo "Starting time:                $(zdump MEC)" >> bwa.log
	echo "Sample file:                  ${sfile}" >> bwa.log
	echo "Reference used for BWA MEM:   ${currentdir}/${ref}" >> bwa.log
	echo "Path to sample:               ${samplepath}/${sample}${ext1}" >> bwa.log
	if [ $mode = 'paired' ]
	then 
	echo "Path to sample [pair]:        ${samplepath}/${sample}${ext2}" >> bwa.log
	fi
	echo "-T bwa mem quality parameter: $qual" >> bwa.log
	echo "-Q mapping quality parameter: $Q" >> bwa.log
	echo "File extension:               $ext" >> bwa.log
	echo "Mode:                         $mode" >> bwa.log
	
	
   	## MAPPING: BWA MEM
	# some bwa mem options
	# t: threads
	# M: flag shorter split hits as secondary (important for piccard tools)
	# w: band width, essentially gaps longer than this are not found (default 100)
	# T: dont output with score lower than T (default 30)
	# B: mismatch penalty (default 4)
	# E: gap extension penalty (default 1)
	
	if [ $mode = 'paired' ]
	then
	
		bwa mem -t 5 -M -T $qual -R "@RG\tID:$sample\tSM:$sample" $ref ${sample}${ext1} ${sample}${ext2} > ${sample}.bwa-mem.sam 2> ${sample}.bwa-mem.log
	
	elif [ $mode = 'unpaired' ]
	then
		
		bwa mem -t 5 -M -T $qual -R "@RG\tID:$sample\tSM:$sample" $ref ${sample}${ext1} > ${sample}.bwa-mem.sam 2> ${sample}.bwa-mem.log
			
	fi
		
	## compress: .sam -> .bam
	#samtools view -b -S ${sample}.bwa-mem.sam -o ${sample}.bwa-mem.bam
        
	## extract mapped reads from raw .sam
	nmapped=$(basename ${sample}.bwa-mem.sam .sam)
	samtools view -S -h  -F4 ${sample}.bwa-mem.sam > ${nmapped}.mapped.sam
	sleep 5
        
	## compress and sort: mapped.sam -> mapped.sorted.bam
	samtools view -b -S ${sample}.bwa-mem.mapped.sam  | samtools sort -o - -  > ${sample}.bwa-mem.mapped.sorted.bam
        
	## index mapped.sorted.bam
	samtools index  ${sample}.bwa-mem.mapped.sorted.bam
	
	## get flagstats
	samtools flagstat ${sample}.bwa-mem.mapped.sorted.bam > ${sample}.bwa-mem.mapped.flagstats &
	 
	## subset .mapped.sam with mapping quality > $Q
	sambase=${sample}.bwa-mem.mapped
	grep "^@" ${sambase}.sam > ${sambase}.Q${Q}.sam
	grep -v "^@"  ${sambase}.sam  |awk -v Q=$Q '{if($5 >= Q ) print $0}'  >> ${sambase}.Q${Q}.sam
	
	## compress and sort: mapped.Q${Q}.sam -> mapped.Q${Q}.sorted.bam 
	newsambase=${sambase}.Q${Q}
	samtools view -b -S ${newsambase}.sam | samtools sort -o - - > ${newsambase}.sorted.bam
        
	## index mapped.Q${Q}.sorted.bam
	samtools index  ${newsambase}.sorted.bam
	
	## get flagstats 
	samtools flagstat ${newsambase}.sorted.bam > ${newsambase}.flagstats &
	
	## get samstats (produces .html and samstat.err file)
	samstat  ${newsambase}.sorted.bam > ${newsambase}.samstat.log &
	
	## Delete .sam files
	/bin/rm -f ${sample}.bwa-mem.sam         
	/bin/rm -f ${sample}.bwa-mem.mapped.sam  
	/bin/rm -f ${sample}.bwa-mem.mapped.Q${Q}.sam 
	#/bin/rm -f ${sample}.bwa-mem.bam
        
	## Sample Finish Time
	echo "Finish time:                  $(zdump MEC)" >> bwa.log
	echo " " >> bwa.log
	
	cd ..
}
    

## Export variables (they are not available inside the doMapping() function otherwise)
export sfile=$sfile
export ref=$ref
export refname=$refname
export qual=$qual
export Q=$Q
export ext=$ext
export ext1=$ext1
export ext2=$ext2
export mode=$mode
export samplepath=$samplepath
export currentdir=$currentdir
export -f doMapping


## DO MAPPING
echo
echo "Start mapping..."
echo
cat ${sfile} | /usr/local/bin/parallel -j $threads doMapping


## Finish
echo 
echo "All samples processed."
echo

