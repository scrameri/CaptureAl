#!/bin/bash

## Usage: select.best.contigs.parallel.sh -s <sample file> -l <locus file> -d <directory with blast results> -c <contigs fasta> -r <reference.fasta> -m <alignment model> -g <FLAG: if given, will not write query and target ranges to .exonerate output> -e <evalue> -p <min perc query in aln> -a <min aln length> -n <min contig length> -x <max contig length> -i <min perc id> -t <number of threads>

## Value: for each sample and locus, gives ${sample}.${locus}.allScore, ${sample}.${locus}.bestScore and ${sample}.${locus}.bestScore.fasta with the best-matching contig (if exonerate could do an alignment)

## Needs:
# exonerate
# awk
# perl
# filter.fasta.by.length.R

## Define arguments
getranges='true'
while getopts s:l:d:c:r:m:e:p:a:n:x:i:t:g opts
do
        case "${opts}"
		in
				s) sfile=${OPTARG};;
				l) lfile=${OPTARG};;
				d) blastdir=${OPTARG};;
				c) contigs=${OPTARG};;
				r) ref=${OPTARG};;
				m) model=${OPTARG};;
				g) getranges='false';;
				e) evalue=${OPTARG};;
				p) minpropaln=${OPTARG};;
				a) minalnlen=${OPTARG};;
				n) minlen=${OPTARG};;
				x) maxlen=${OPTARG};;
				i) minpercid=${OPTARG};;
				t) threads=${OPTARG};;
		esac
done

## Check parameters and set defaults
if [ ! $sfile ] ; then echo "sample file not specified (-s option), stopping" ; exit 0 ; fi
if [ ! $lfile ] ; then echo "locus file not specified (-l option), stopping" ; exit 0 ; fi
if [ ! $blastdir ] ; then echo "directory with blast results not specified (-d option), setting to 'blast.contigs'" ; blastdir=blast.contigs ; fi
if [ ! $contigs ] ; then echo "name of contigs fasta in sample subdirs not specified (-c option), stopping" ; exit 0 ; fi
if [ ! $ref ] ; then echo "reference .fasta not specified (-r option), stopping" ; exit 0 ; fi

if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping" ; exit 0 ; fi
if [ ! -f $lfile ] ; then echo "locus file <$lfile> not found, stopping" ; exit 0 ; fi
if [ ! -d $blastdir ] ; then echo "directory with blast results <$blastdir> not found, stopping" ; exit 0 ; fi
if [ ! -f $ref ] ; then echo "reference .fasta <$ref> not found, stopping" ; exit 0 ; fi

if [ $getranges == 'true' ] ; then echo "will write query and target ranges to .exonerate output (takes some more time)" ; fi

if [ ! $evalue ] ; then echo "evalue threshold (-e option) not set, looking for .blast files generated using evalue threshold of 1e-04" ; evalue=1e-04 ; fi

if [ ! $minpropaln ] ; then echo "minimum minimum percentage of query.length in alignment length (-p option) not set, setting to 0" ; minpropaln=0 ; fi
if [ ! $minalnlen ] ; then echo "minimum alignment length (-a option) not set, setting to 80" ; minalnlen=80 ; fi
if [ ! $minlen ] ; then echo "minimum contig length (-n option) not set, setting to 100" ; maxlen=100 ; fi
if [ ! $maxlen ] ; then echo "maximum contig length (-x option) not set, setting to 10000" ; maxlen=10000 ; fi
if [ ! $minpercid ] ; then echo "minimum percent identity (-i option) not set, setting to 80" ; minpercid=80 ; fi

if [ ! $threads ] ; then echo "number of threads (-t option) not set, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30" ; exit 0 ; fi

## Available gapped alignment options (from exonerate man page) (-m option)
# - affine:global
#   This performs gapped global alignment, similar to the Needleman-Wunsch algorithm, except with affine gaps. Global alignment requires that both the sequences in their entirety are included in the alignment.
# - affine:bestfit
#   This performs a best fit or best location alignment of the query onto the target sequence. The entire query sequence will be included in the alignment, but only the best location for its alignment on the target sequence.
# - affine:local
#   This is local alignment with affine gaps, similar to the Smith-Waterman-Gotoh algorithm. A general-purpose alignment algorithm. As this is local alignment, any subsequence of the query and target sequence may appear in the alignment.
if [ ! $model ] ; then echo "gapped alignment model not specified (-m option), setting to 'affine:local'" ; model='affine:local' ; fi
if [ $model != 'affine:local' ] && [ $model != 'affine:global' ] && [ $model != 'affine:bestfit' ] ; then echo "gapped alignment model <$model> not supported (choose between 'affine:local', 'affine:global' or 'affine:bestfit'" ; exit 0 ; fi

## Export parameters
export sfile=$sfile
export lfile=$lfile
export blastdir=$blastdir

export ref=$ref
export refbase=$(basename $ref .fasta)

export model=$model

export evalue=$evalue
export minpropaln=$minpropaln
export minalnlen=$minalnlen
export minlen=$minlen
export maxlen=$maxlen
export minpercid=$minpercid
nind=$(wc -l < ${sfile})
nloc=$(wc -l < ${lfile})


## Additional arguments
export dir=$(pwd)
export outdir="best.contigs.${nind}.${nloc}"  # name of output directory
export contigbase=$(basename $contigs .fasta) # contig file in $blastdir/$sample is $contigbase.fasta
export resultempty="--------------------"     # encoding for missing sequences (no contig or no alignment at a locus)
export getranges=$getranges                   # if 'true' will write .exonerate with <query> <target> <score> <query range> <target range> outputs

## Extractor function (extracts values after a string in exonerate results)
do_extract()
{
file=$1
string="$2 "
pstring=" $2"
sstring="^$(echo $pstring | awk '/,/(gsub(/ /,""))')"
grep "$string" ${file} | awk '/,/(gsub(/ /,""))' | awk -v str=$sstring '/,/(gsub(str,""))'
}
export -f do_extract

## Check that .blast files for each sample are available
nfiles=$(ls -1 ${dir}/${blastdir}/*/${refbase}.on.*.${contigbase}.${evalue}.${minpropaln}.${minalnlen}.${minlen}.${maxlen}.${minpercid}.passedhits.txt | wc -l) 
if [ $nind -gt $nfiles ] ; then echo "number of individuals in <${sfile}> is ${nind}, but number of <.blast> files is only ${nfiles}, stopping." ; exit 0 ; fi

## Be verbose
echo
echo "looking for BLAST hits in ${blastdir}/{sample}/${refbase}.on.{sample}.${contigbase}.${evalue}.${minpropaln}.${minalnlen}.${minlen}.${maxlen}.${minpercid}.passedhits.txt"
echo 

## Generate log file
echo ""                                                                  > best.contigs.${nind}.${nloc}.log
echo "=== select.best.contigs LOG ==="                                   >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "sample file:                                         $sfile"       >> best.contigs.${nind}.${nloc}.log
echo "locus file:                                          $lfile"       >> best.contigs.${nind}.${nloc}.log
echo "blast directory:                                     $blastdir"    >> best.contigs.${nind}.${nloc}.log
echo "selects best matching contigs from this file:        $contigs"     >> best.contigs.${nind}.${nloc}.log
echo "matches contigs to this reference:                   $ref"         >> best.contigs.${nind}.${nloc}.log
echo "exonerate gapped alignment model:                    $model"       >> best.contigs.${nind}.${nloc}.log
echo "write query and target ranges to .exonerate output:  $getranges"   >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "looking for BLAST hits based on evalue threshold:    $evalue"      >> best.contigs.${nind}.${nloc}.log
echo "minimum percent of query length in BLAST alignment:  $minpropaln"  >> best.contigs.${nind}.${nloc}.log
echo "minimum BLAST alignment length:                      $minalnlen"   >> best.contigs.${nind}.${nloc}.log
echo "minimum contig length:                               $minlen"      >> best.contigs.${nind}.${nloc}.log
echo "maximum contig length:                               $maxlen"      >> best.contigs.${nind}.${nloc}.log
echo "minimum percent identity in BLAST alignment:         $minpercid"   >> best.contigs.${nind}.${nloc}.log
echo "number of threads:                                   $threads"     >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "number of loci to work on:                           $nloc"        >> best.contigs.${nind}.${nloc}.log
echo "number of samples to work on:                        $nind"        >> best.contigs.${nind}.${nloc}.log
echo "number of BLAST files:                               $nfiles"      >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "Starting time:                                       $(zdump MEC)" >> best.contigs.${nind}.${nloc}.log


## Create output directory
if [ -d ${outdir} ] 
then
	mv ${outdir} ${outdir}.bak ; echo ; echo "${outdir} exists, moving it to ${outdir}.bak" ; echo
fi

mkdir ${outdir}

## Extract reference sequences
cd ${outdir}
echo ; echo "extracting reference sequences..."

# Index and extract reference sequences
#mkdir refseqs
#cd refseqs
#
#fastaindex ${dir}/${ref} ${ref}.idx
#
#for locus in $(cat ${dir}/${lfile})
#do
#	fastafetch -f ${dir}/${ref} -i ${ref}.idx -q ${locus} > ${locus}.fasta
#done	
extract.all.fasta.seqs.R ${dir}/${ref} refseqs > /dev/null

echo "done" ; echo


## Define function to select best contig using exonerate
get_best_contig()
{
	sample=$1
	echo "${sample}"
	#sampleBase=$(basename ${sample} ${suffix})
	
	mkdir ${sample}
	cd ${sample}
		
	# index contigs
	allcontigs="${dir}/${blastdir}/${sample}/${contigbase}.fasta"
	if [ ! -f $allcontigs ] ; then echo "$allcontigs not found" >> ${dir}/best.contigs.${nind}.${nloc}.err ; fi
	# allcontigs file looks like: blast.dipspades.contigs/RIR2458_S13_L001/consensus_contigs.fasta
	fastaindex ${allcontigs} ${sample}.allcontigs.idx
	
	# path to .blast file
	#blastfile="${dir}/${blastdir}/${sample}/${refbase}.on.${sample}.${contigbase}.${evalue}.blast" # all hits (contig name in column 3)
	blastfile="${dir}/${blastdir}/${sample}/${refbase}.on.${sample}.${contigbase}.${evalue}.${minpropaln}.${minalnlen}.${minlen}.${maxlen}.${minpercid}.passedhits.txt" # passed hits (contig name in column 5)
	if [ ! -f $blastfile ] ; then echo "$blastfile not found" >> ${dir}/best.contigs.${nind}.${nloc}.err ; fi	
		
	for locus in $(cat ${dir}/${lfile})
	do 
			
		# fetch contigs with .blast hit to locus
		#grep -w "^$locus" ${blastfile} | cut -f3 | sort | uniq > ${sample}.${locus}.contignames
		grep -w "$locus" ${blastfile} | cut -f5 | sort | uniq > ${sample}.${locus}.contignames
		fastafetch -f ${allcontigs} -i ${sample}.allcontigs.idx -Fq ${sample}.${locus}.contignames > ${sample}.${locus}.contigs.fasta
		
		# discard contigs shorter than $minlen or longer than $maxlen (if input is not already filtered)
		#filter.fasta.by.length.R ${sample}.${locus}.contigs.fasta $minlen $maxlen
		#mv ${sample}.${locus}.contigs_l${minlen}-${maxlen}.fasta ${sample}.${locus}.contigs.fasta
		
		# export a.fa (contigs) and b.fa (reference)
		export contigs="${sample}.${locus}.contigs.fasta"
		export reference="${dir}/${outdir}/refseqs/${locus}.fasta"
		
		# run exonerate
		if [ $getranges == 'true' ]
		then
			# slower, but also outputs query, target, raw score, query range and target range
			exonerate --model ${model} --bestn 1 --exhaustive yes ${contigs} ${reference} --showcigar FALSE --showvulgar FALSE --showsugar TRUE --showalignment FALSE 2> /dev/null | \
				grep '^sugar:' | sed -e 's/^sugar: //g' | sort -k9nr > ${sample}.${locus}.exonerate
			#looks like: 
			#5983_contig_716_length 596 371 - consFabaceae_NA_LG_Scaffold084976_1_150_ID_2979 1018 1239 + 104
			#7953_contig_653_length 610 173 - consFabaceae_NA_LG_Scaffold084976_1_150_ID_2979 1294 1710 + 123
			#corresponds to: <1: query> <2: qstart> <3: qend> <4: -> <5: target> <6: tstart> <7: tend> <8: +> <9: raw score>
			result=$(cut -f1,9 -d' ' ${sample}.${locus}.exonerate | xargs)
			#looks like:
			#5983_contig_716_length 104 7953_contig_653_length 123
			
			#exonerate --model ${model} --bestn 1 --exhaustive yes ${contigs} ${reference} 1> ${sample}.${locus}.exonerate 2>/dev/null
			#query=$(do_extract ${sample}.${locus}.exonerate "Query:" | sed -e "s/\[revcomp\]$//g")
			#target=$(do_extract ${sample}.${locus}.exonerate "Target:")
			#score=$(do_extract ${sample}.${locus}.exonerate "Raw score:")
			#qrange=$(do_extract ${sample}.${locus}.exonerate "Query range:")
			#trange=$(do_extract ${sample}.${locus}.exonerate "Target range:")
			#paste <(printf %.0s"$sample\n" $(seq $(echo $query | wc -w))) <(printf %s "$query") <(printf %s "$target") <(printf %s "$score") <(printf %s "$qrange") <(printf %s "$trange") | sort -k4nr > ${sample}.${locus}.exonerate
			#result=$(cut -f2,4 ${sample}.${locus}.exonerate | xargs)
			
		else
			# faster, but only outputs query and raw score
			result=`perl -w -e 'system("exonerate --model $ENV{model} --bestn 1 --exhaustive yes  $ENV{contigs}  $ENV{reference}  2>/dev/null |grep -P \"Query:|Raw score\" |perl -p -e \"s/Query\:|Raw score\:|revcomp//g\" | tr -d \"[]\"  ")' `
			# result looks like: 5983_contig_716_length 104 7953_contig_653_length 123
		fi
		
		# write main output files
		if [ -z "${result}" ]  # if $result is empty, can happen if exnoerate is not able to produce useful alignment
		then 
			echo "0" > ${sample}.${locus}.allScore
			echo "0" > ${sample}.${locus}.bestScore
			
			echo ">${sample}" > ${sample}.${locus}.bestScore.fasta
			echo ${resultempty} >> ${sample}.${locus}.bestScore.fasta
		else 
			read bestscore bestname <<< `echo $result | perl -p -e 's/\s+/\n/g' | \
				perl -w -e 'my $item=0; my $bestscore=0; while ($l = <>){chomp $l; $item++; if($item %2 ==1){$name=$l};  if($item %2 ==0){if ($l > $bestscore) {$bestscore = $l ; $bestname = $name }}} print "$bestscore  $bestname";'`
			echo ${result} > ${sample}.${locus}.allScore
			bestnamescore=$(echo $bestname $bestscore)
			echo "$bestnamescore" > ${sample}.${locus}.bestScore
			
			echo ">${sample}" > ${sample}.${locus}.bestScore.fasta
			fastafetch -f ${allcontigs} -i ${sample}.allcontigs.idx -q ${bestname} | \
				awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | grep -v '^>' >> ${sample}.${locus}.bestScore.fasta
		fi
		
		# clean up
		rm -f ${sample}.${locus}.contignames
		rm -f ${sample}.${locus}.contigs.fasta
	
	done
	
	rm -f ${sample}.allcontigs.idx
	cd ${dir}/${outdir}
}


## Execute function
export -f get_best_contig
cd ${dir}/${outdir}
cat ${dir}/${sfile} | parallel -j $threads get_best_contig
cd ${dir}


## Finish
echo 
echo "Finish time:                                         $(zdump MEC)" >> best.contigs.${nind}.${nloc}.log
echo "All samples processed."
echo
