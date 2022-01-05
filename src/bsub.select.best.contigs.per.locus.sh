#!/bin/bash

#BSUB -J "SEL[1-15]%5"
#BSUB -R "rusage[mem=250]"          # 200 to 250 are sufficient 
#BSUB -n 1
#BSUB -W 4:00

## Usage: bsub < bsub.select.best.contigs.per.locus.sh # from spades-2396/run directory

## Needs: exonerate, extract.all.fasta.seqs.R

# load modules
module load gdc exonerate/2.4.0 r/3.1.2

# arguments
in=$(pwd)
wd=$(pwd)
sfile="samples.txt"
reg="testloci.txt"
ref="/cluster/home/crameris/home/Dalbergia/uce/references/consDalbergia_4c_2396.fasta"
out="/cluster/scratch/crameris/combined/best.contigs.$(wc -l < ${sfile}).$(wc -l < ${reg})"

# additional arguments
prefix="extracted_reads_"
suffix=".spades"
contigspath="SAMPLE/${prefix}SAMPLE.LOCUS${suffix}/contigs.fasta" # relative to ${in}
getranges="true"
model="affine:local"
exhaustive="FALSE" # "FALSE" or "yes" (SLOW!) 
resultempty="--------------------"
#threads=2

# check arguments
if [ ! -f ${sfile} ] ; then echo "ERROR: sample file <${sfile}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${reg} ] ; then echo "ERROR: locus file <${reg}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${ref} ] ; then echo "ERROR: reference .fasta file <${ref}> not found, stopping" ; exit 0 ; fi
if [ ! -d ${in} ] ; then echo "ERROR: input directory (reads) <${in}> not found, stopping" ; exit 0 ; fi

# exonerate options
# available gapped alignment options (from exonerate man page) (-m option)
# - affine:global
#   This performs gapped global alignment, similar to the Needleman-Wunsch algorithm, except with affine gaps. Global alignment requires that both the sequences in their entirety are included in the alignment.
# - affine:bestfit
#   This performs a best fit or best location alignment of the query onto the target sequence. The entire query sequence will be included in the alignment, but only the best location for its alignment on the target sequence.
# - affine:local
#   This is local alignment with affine gaps, similar to the Smith-Waterman-Gotoh algorithm. A general-purpose alignment algorithm. As this is local alignment, any subsequence of the query and target sequence may appear in the alignment.
if [ ! ${model} ] ; then echo "exonerate gapped alignment model not specified (--model option), setting to 'affine:local'" ; model='affine:local' ; fi
if [ ${model} != 'affine:local' ] && [ ${model} != 'affine:global' ] && [ ${model} != 'affine:bestfit' ] ; then echo "gapped alignment model <$model> not supported (choose between 'affine:local', 'affine:global' or 'affine:bestfit'" ; exit 0 ; fi

if [ ! ${exhaustive} ] ; then echo "exonerate --exhaustive option not specified, setting to 'FALSE'" ; exhaustive='FALSE' ; fi
if [ ${exhaustive} != 'FALSE' ] && [ ${model} != 'yes' ] ; then echo "exonerate --exhaustive option not supported (choose between 'FALSE' or 'yes'" ; exit 0 ; fi


# get number of samples and regions
nind=$(wc -l < ${sfile})
nloc=$(wc -l < ${reg})

# create output directory
if [ ! -d $(dirname ${out}) ] ; then mkdir $(dirname ${out}) ; fi
if [ ! -d ${out} ] ; then mkdir ${out} ; fi

if [ ! -f ${out}/$(basename ${sfile}) ] ; then cp ${sfile} ${out} ; fi
if [ ! -f ${out}/$(basename ${reg}) ] ; then cp ${reg} ${out} ; fi

# get job index
IDX=$LSB_JOBINDEX
name=$(sed -n ${IDX}p < ${sfile})

# extract reference sequences
if [ -d ${out}/refseqs ] ; then /bin/rm -rf ${out}/refseqs ; fi
extract.all.fasta.seqs.R ${ref} ${out}/refseqs > /dev/null

# create sample subdirectory)
cd ${out}
if [ ! -d ${name} ] ; then mkdir ${name} ; fi
cpath=$(echo ${contigspath} | sed -e "s/SAMPLE/${name}/g")

# extract assemblies from input archive ${inarch}
inarch=$(echo $(dirname $(dirname ${in}/${contigspath})).tar.gz | sed -e "s/SAMPLE/${name}/g")
tar -xf ${inarch}
if [ -f ${name}/assembly.log ] ; then /bin/rm -f ${name}/assembly.log ; fi
if [ -f ${name}/${name}.cmds.txt ] ; then /bin/rm -f ${name}/${name}.cmds.txt ; fi

# generate log file
logfile="${name}/best.contigs.${nind}.${nloc}.log"

echo "=================================================================" > ${logfile}
echo "========================== ASSEMBLY LOG =========================" >> ${logfile}
echo "=================================================================" >> ${logfile}
echo " " >> ${logfile}
echo "Starting time:                    $(zdump CET)" >> ${logfile}
echo "Working directory:                ${wd}" >> ${logfile}
echo "Input directory:                  ${in}" >> ${logfile}
echo "Output directory:                 ${out}" >> ${logfile}
echo "Sample file:                      ${sfile} (${nind} samples)" >> ${logfile}
echo "Locus file:                       ${reg} (${nloc} regions)" >> ${logfile}
echo "Input archive:                    ${inarch}" >> ${logfile}
echo "Contigs file (from sample dir):   ${cpath}" >> ${logfile}
echo "Reference .fasta file:            ${ref}" >> ${logfile}
echo "Exonerate gapped alignment model: ${model}" >> ${logfile}
echo "Write query and target ranges:    ${getranges}" >> ${logfile}
echo "Exhaustive alignment:             ${exhaustive}" >> ${logfile}
#echo "Number of threads used:       ${threads}" >> ${logfile}	

# loop through regions
for region in $(cat ${out}/$(basename ${reg}))
do
	
	cd ${name}
	
	# export a.fa (contigs) and b.fa (reference)
	export cpath=$(echo "${contigspath}" | sed -e "s/SAMPLE/${name}/g" | sed -e "s/LOCUS/${region}/g")
	export reference="${out}/refseqs/${region}.fasta"

	# peek into .tar
	#tar -tvf ${prefix}${name}.${region}${suffix}.tar.gz

	# extract assembled contigs
	#tar -xf ${prefix}${name}.${region}${suffix}.tar.gz --strip-components=1 ${name}/${prefix}${name}.${region}${suffix}/spades.log
	tar -xf ${prefix}${name}.${region}${suffix}.tar.gz --strip-components=1 ${cpath} 2> /dev/null

	contigs=$(echo $(basename $(dirname $contigs))/$(basename $contigs))	
	
	# index contigs
	idx="${prefix}${name}.${region}${suffix}/contigs.idx"
	fastaindex ${contigs} ${idx} 2> /dev/null
	
	# run exonerate
	# NOTE: in exonerate version 2.4 and using spades with --isolate, exonerate is not able to produce --exhaustive yes alignments anymore???
	
	if [ $getranges == 'true' ]
	then
		# slower, but also outputs query, target, raw score, query range and target range
		exonerate --model ${model} --bestn 1 --exhaustive FALSE ${contigs} ${reference} --showcigar FALSE --showvulgar FALSE --showsugar TRUE --showalignment FALSE 2> /dev/null | \
			grep '^sugar:' | sed -e 's/^sugar: //g' | sort -k9nr > ${name}.${region}.exonerate
		#looks like: 
		#5983_contig_716_length 596 371 - consFabaceae_NA_LG_Scaffold084976_1_150_ID_2979 1018 1239 + 104
		#7953_contig_653_length 610 173 - consFabaceae_NA_LG_Scaffold084976_1_150_ID_2979 1294 1710 + 123
		#corresponds to: <1: query> <2: qstart> <3: qend> <4: -> <5: target> <6: tstart> <7: tend> <8: +> <9: raw score>
		result=$(cut -f1,9 -d' ' ${name}.${region}.exonerate | xargs)
		#looks like:
		#5983_contig_716_length 104 7953_contig_653_length 123
		
		#exonerate --model ${model} --bestn 1 --exhaustive yes ${contigs} ${reference} 1> ${name}.${region}.exonerate 2>/dev/null
		#query=$(do_extract ${name}.${region}.exonerate "Query:" | sed -e "s/\[revcomp\]$//g")
		#target=$(do_extract ${name}.${region}.exonerate "Target:")
		#score=$(do_extract ${name}.${region}.exonerate "Raw score:")
		#qrange=$(do_extract ${name}.${region}.exonerate "Query range:")
		#trange=$(do_extract ${name}.${region}.exonerate "Target range:")
		#paste <(printf %.0s"$name\n" $(seq $(echo $query | wc -w))) <(printf %s "$query") <(printf %s "$target") <(printf %s "$score") <(printf %s "$qrange") <(printf %s "$trange") | sort -k4nr > ${name}.${region}.exonerate
		#result=$(cut -f2,4 ${name}.${region}.exonerate | xargs)
		
	else
		# faster, but only outputs query and raw score
		result=$(exonerate --model ${model} --bestn 1 --exhaustive FALSE --showalignment yes ${contigs}  ${reference} 2>/dev/null | \
			grep -P 'Query:|Raw score' |Rscript -e 'x <- sapply(strsplit(readLines("stdin"), split = ": "), "[", 2) ; y <- data.frame(c=x[seq(1,length(x),by=2)],s=x[seq(2,length(x),by=2)],stringsAsFactors=F) ; z <- y[order(as.numeric(y$s),decreasing=T),] ; cat(paste(as.character(t(as.matrix(z))), collapse = " "))' )
		# result looks like: 5983_contig_716_length 104 7953_contig_653_length 123
	fi
	
	# write main output files
	if [ -z "${result}" ]  # if $result is empty, can happen if exnoerate is not able to produce useful alignment
	then 
		echo "0" > ${name}.${region}.allScore
		#echo "0" > ${name}.${region}.bestScore
		
		echo ">${name}" > ${name}.${region}.bestScore.fasta
		echo ${resultempty} >> ${name}.${region}.bestScore.fasta
	else 
		echo ${result} > ${name}.${region}.allScore
		
		# old way (writes unnecessary ${name}.${region}.bestScore file)
		#read bestscore bestname <<< $(echo $result | perl -p -e 's/\s+/\n/g' | perl -w -e 'my $item=0; my $bestscore=0; while ($l = <>){chomp $l; $item++; if($item %2 ==1){$name=$l};  if($item %2 ==0){if ($l > $bestscore) {$bestscore = $l ; $bestname = $name }}} print "$bestscore  $bestname";'`
		#bestnamescore=$(echo $bestname $bestscore)
		#echo "$bestnamescore" > ${name}.${region}.bestScore
		
		# new way (assumes ${result} is already sorted for decreasing score)
		bestname=$(echo ${result} | cut -f1 -d' ')
		echo ">${name}" > ${name}.${region}.bestScore.fasta
		fastafetch -f ${contigs} -i ${idx} -q ${bestname} | \
			awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | grep -v '^>' >> ${name}.${region}.bestScore.fasta 
	fi
		
	# clean up
	/bin/rm -rf ${prefix}${name}.${region}${suffix}
	/bin/rm -f ${prefix}${name}.${region}${suffix}.tar.gz
	cd ${out}
done	

# clean up
/bin/rm -f ${name}/*${suffix}.tar.gz

# finish
echo "Finish time:                      $(zdump CET)" >> ${logfile}

# tar.gz sample subdirectory with .allScore / .bestScore.fasta (.exonerate) results
tar -czf ${name}.tar.gz ${name}
/bin/rm -rf ${name}
