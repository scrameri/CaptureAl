#!/bin/bash

## Usage: select.best.contigs.per.locus.sh -s <samples.txt> -l <loci.txt> -r <reference.fasta> -d <input directory> -c <path from input directory to contigs, will replace 'SAMPLE' and 'LOCUS' with the respective strings> -g <FLAG: if given, will not write query and target ranges to .exonerate output> -m <alignment model> -t <number of threads>

## Needs: exonerate, extract.all.fasta.seqs.R

## Define arguments
getranges='true'
while getopts s:l:r:d:c:g:m:t: opts
do
        case "${opts}"
        in
        	s) sfile=${OPTARG};;
        	l) lfile=${OPTARG};;
        	r) ref=${OPTARG};;
			d) indir=${OPTARG};;
			c) contigspath=${OPTARG};;
			g) getranges='false';;
			m) model=${OPTARG};;
			t) threads=${OPTARG};;
    	esac
done

## Check input
if [ ! $sfile ] ; then echo "sample file not provided (-s option), stopping." ; exit 0 ; fi
if [ ! $lfile ] ; then echo "locus file not provided (-l option), stopping." ; exit 0 ; fi
if [ ! $ref ] ; then echo "reference sequences not provided (-r option), stopping." ; exit 0 ; fi
if [ ! $indir ] ; then echo "input directory not provided (-d option), stopping." ; exit 0 ; fi
if [ ! $contigspath ] ; then echo "path from input directory to contigs not provided (-c option), setting to <SAMPLE.dipspades/extracted_reads_SAMPLE.fastq.LOCUS.ids.spades/dipspades/consensus_contigs.fasta>." ; contigspath="SAMPLE.dipspades/extracted_reads_SAMPLE.fastq.LOCUS.ids.spades/dipspades/consensus_contigs.fasta" ; fi
if [ $getranges == 'true' ] ; then echo "will write query and target ranges to .exonerate output (takes some more time)" ; fi
if [ ! -f $sfile ] ; then echo "sample file <$sfile> not found, stopping." ; exit 0 ; fi
if [ ! -f $lfile ] ; then echo "locus file <$lfile> not found, stopping." ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "input directory <$indir> not found, stopping." ; exit 0 ; fi
if [ ! -f $ref ] ; then echo "reference sequences <$ref> not found, stopping." ; exit 0 ; fi
if [ ! $threads ] ; then threads=15 ; fi

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
export indir=${indir}

export ref=$ref
export refbase=$(basename $ref .fasta)

export contigspath=${contigspath}
export model=$model
nind=$(wc -l < ${sfile})
nloc=$(wc -l < ${lfile})

## Additional arguments
export dir=$(pwd)
export outdir="best.contigs.${nind}.${nloc}"  # name of output directory
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

## Create output directory
if [ -d ${outdir} ] 
then
	mv ${outdir} ${outdir}.bak ; echo ; echo "${outdir} exists, moving it to ${outdir}.bak" ; echo
fi

mkdir ${outdir}

## Get sample names
cp ${sfile} ${outdir}/samples.txt

## Get locus names
cp ${lfile} ${outdir}/loci.txt

## Extract reference sequences
cd ${outdir}
echo ; echo "extracting reference sequences..."
extract.all.fasta.seqs.R ${dir}/${ref} refseqs > /dev/null
echo "done" ; echo ; cd ${dir}

## Generate log file
echo ""                                                                  > best.contigs.${nind}.${nloc}.log
echo "=== select.best.contigs LOG ==="                                   >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "sample file:                                         $sfile"       >> best.contigs.${nind}.${nloc}.log
echo "locus file:                                          $lfile"       >> best.contigs.${nind}.${nloc}.log
echo "selects best matching contigs from this file:        $contigspath" >> best.contigs.${nind}.${nloc}.log
echo "matches contigs to this reference:                   $ref"         >> best.contigs.${nind}.${nloc}.log
echo "exonerate gapped alignment model:                    $model"       >> best.contigs.${nind}.${nloc}.log
echo "write query and target ranges to .exonerate output:  $getranges"   >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "number of threads:                                   $threads"     >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "number of loci to work on:                           $nloc"        >> best.contigs.${nind}.${nloc}.log
echo "number of samples to work on:                        $nind"        >> best.contigs.${nind}.${nloc}.log
echo ""                                                                  >> best.contigs.${nind}.${nloc}.log
echo "Starting time:                                       $(zdump MEC)" >> best.contigs.${nind}.${nloc}.log

## Define function
get_best_contig()
{
	sample=$1
	echo "$sample"
	
	mkdir ${sample}
	cd ${sample}	

	for locus in $(cat ../loci.txt)
	do
		
		# export a.fa (contigs) and b.fa (reference)
		export contigs=$(echo "${dir}/${indir}/${contigspath}" | sed -e "s/SAMPLE/$sample/g" | sed -e "s/LOCUS/$locus/g")
		export reference="${dir}/${outdir}/refseqs/${locus}.fasta"
		
		# index contigs
		fastaindex ${contigs} ${sample}.contigs.idx 2> /dev/null
		
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
			read bestscore bestname <<< $(echo $result | perl -p -e 's/\s+/\n/g' | perl -w -e 'my $item=0; my $bestscore=0; while ($l = <>){chomp $l; $item++; if($item %2 ==1){$name=$l};  if($item %2 ==0){if ($l > $bestscore) {$bestscore = $l ; $bestname = $name }}} print "$bestscore  $bestname"')
			echo ${result} > ${sample}.${locus}.allScore
			bestnamescore=$(echo $bestname $bestscore)
			echo "$bestnamescore" > ${sample}.${locus}.bestScore
			
			echo ">${sample}" > ${sample}.${locus}.bestScore.fasta
			fastafetch -f ${contigs} -i ${sample}.contigs.idx -q ${bestname} | \
				awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | grep -v '^>' >> ${sample}.${locus}.bestScore.fasta 
		fi
		
		# clean up
		rm -f ${sample}.contigs.idx
	
	done	
	cd ${dir}/${outdir}
}


## Execute function
export -f get_best_contig
cd ${dir}/${outdir}
cat samples.txt | parallel -j $threads get_best_contig
cd ${dir}


## Finish
echo 
echo "Finish time:                                         $(zdump MEC)" >> best.contigs.${nind}.${nloc}.log
echo "All samples processed."
echo
