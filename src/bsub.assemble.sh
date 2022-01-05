#!/bin/bash

#BSUB -J "ASS[1,5,10]%3"            # [1-45]%15 
#BSUB -R "rusage[mem=250]"          # 200 to 250 are sufficient 
#BSUB -n 1							# leave at 1. 2 threads used for hyperthreading
#BSUB -W 24:00                      # 24:00 needed

## Usage: bsub < bsub.assemble.sh # from extracted.Q${Q}/run directory

# load modules
module load new gcc/6.3.0 gdc python/2.7.11 spades/3.14.1

# arguments
in=$(pwd) # need absolute path
out="/cluster/scratch/crameris/spades-2396/$(basename ${in})"
wd=$(pwd) # working directory
threads=2 # number of threads for hyperthreading
sfile="samples.txt"
reg="loci.txt" # ${in}.loci.txt

# additional arguments
prefix="extracted_reads_" # prefix of forward read .fastq filename
infix1=".trim1." # infix of forward read .fastq filename
infix2=".trim2." # infix of reverse read .fastq filename
suffix="fastq" # suffix of forward read .fastq filename
infixbase="trim1.fastq." # will be removed from $readbase1 to construct $readbase

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${reg} ] ; then echo "ERROR: locus file <${reg}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${in} ] ; then echo "ERROR: input directory (extracted reads) <${in}> not found, stopping" ; exit 0 ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# create output directory
if [ ! -d $(dirname ${out}) ] ; then mkdir $(dirname ${out}) ; fi
if [ ! -d ${out} ] ; then mkdir ${out} ; fi

# copy samples and regions file
cp ${sfile} ${out}
cp ${reg} ${out}

# change to output directory
cd ${out}

# create and clean up sample subdirectory
if [ ! -d ${name} ] ; then mkdir ${name} ; fi
reads="${in}/${name}" # sample subdirectory in extracted reads directory
/bin/rm -rf ${name}/${name}.cmds.txt

# create log file
logfile="${name}/assembly.log"

echo "=================================================================" > ${logfile}
echo "========================== ASSEMBLY LOG =========================" >> ${logfile}
echo "=================================================================" >> ${logfile}
echo " " >> ${logfile}
echo "Starting time:                $(zdump CET)" >> ${logfile}
echo "Working directory:            ${wd}" >> ${logfile}
echo "Input directory:              ${in}" >> ${logfile}
echo "Input read pairs used:        ${reads}.tar.gz" >> ${logfile}
echo "Output directory:             ${out}/${name}" >> ${logfile}
echo "Sample file:                  ${sfile}" >> ${logfile}
echo "Locus file:                   ${reg}" >> ${logfile}
echo "Number of threads used:       ${threads}" >> ${logfile}

# extract reads from .tar.gz file (or copy from extracted directory)
if [ -d ${reads} ]
then
	# copy all .fastq files (will have new timestamp)	
	cp ${reads}/*${suffix} ${name}
elif [ -f ${reads}.tar.gz ]
then
	# extract all .fastq files (will have original timestamp)
	tar -xf ${reads}.tar.gz
	/bin/rm -f ${name}/doExtract.log
else
	echo "${reads} / ${reads}.tar.gz not found, stopping" 
	exit 0
fi

# prepare assembly jobs in loop (over loci)
for region in $(cat $(basename ${reg}))
do
	# example ${fastq1}: ${name}/extracted_reads_ALR2408_L2.trim1.fastq.NW_017997114.1_LG_Scaffold115413_229_328_ID_3906.ids.fastq	
	fastq1="${prefix}${name}${infix1}${suffix}.${region}.ids.${suffix}"
	fastq2="${prefix}${name}${infix2}${suffix}.${region}.ids.${suffix}"
	readbase=$(echo ${fastq1} |sed -e "s/${infixbase}//" |sed -e "s/.ids.${suffix}//")

	# run SPAdes with cov-cutoff auto + careful (with read error correction)
	# If you have high-coverage data for bacterial/viral isolate or multi-cell organism, we highly recommend to use --isolate option
	# to improve run time, use flag --isolate (incompatible with --careful) for high coverage data
	# to omit read error correction, use flag --only-assembler
	# to deactivate auto coverage cutoff, use flag --cov-cutoff off
	# to set a temp dir, use flag --tmp-dir ${TMPDIR}
	# to run dipspades.py, replace spades.py with dipspades executable
	# to use a differernt assembler, change cmd= definition
	cmd="spades.py -t 1 --pe1-1 ${name}/${fastq1} --pe1-2 ${name}/${fastq2} -o ${name}/${readbase}.spades --isolate --cov-cutoff off > /dev/null " #2> ${name}/${readbase}.spades.err"

	echo ${cmd}

done >> ${name}/${name}.cmds.txt

# run assembly jobs in parallel (over loci)
doASS() {
	cmd=$1

	odir=$(echo ${cmd} | cut -f9 -d' ')
	fastq1=$(echo $cmd | cut -f5 -d' ')
	fastq2=$(echo $cmd | cut -f7 -d' ')

	# run assembly command
	if [[ -d ${odir} ]] ; then /bin/rm -rf ${odir} ; fi
	if [[ -f ${odir}.tar.gz ]] ; then /bin/rm -f ${odir}.tar.gz ; fi
	eval "${cmd}"
	
	# clean up output (save disk space, depends on assembler)
	/bin/rm -rf ${fastq1}
	/bin/rm -rf ${fastq2}
	/bin/rm -rf ${odir}/tmp
	/bin/rm -rf ${odir}/K*
	
	# clean up tempfiles written to scratch (by SPAdes assembler, no idea why this happens)
	/bin/rm -rf $(dirname $(dirname ${out}))/core_spades-*
	
	# extract assembled contigs and scaffolds (depends on assembler)
	#if [ -f ${odir}/contigs.fasta ] ; then cp ${odir}/contigs.fasta $(dirname ${odir}) ; fi
	#if [ -f ${odir}/scaffolds.fasta ] ; then cp ${odir}/scaffolds.fasta $(dirname ${odir}) ; fi
	#cp ${odir}/spades.log $(dirname ${odir})

	# tar assembly results
	tar -czf ${odir}.tar.gz ${odir}
	/bin/rm -rf ${odir}
}
export out=${out}
export -f doASS
cat ${name}/${name}.cmds.txt | parallel -j ${threads} doASS

# sample Finish Time
echo "Assembler commands used:      ${name}/${name}.cmds.txt" >> ${logfile}
echo "Finish time:                  $(zdump CET)" >> ${logfile}
echo " " >> ${logfile}

# tar.gz sample subdirectory with (compressed) .spades results
tar -czf ${name}.tar.gz ${name}
/bin/rm -rf ${name}

