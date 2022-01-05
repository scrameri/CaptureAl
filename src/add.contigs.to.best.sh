#!/bin/bash

## Usage: add.contigs.to.best.sh -r <folder with sequences of 1 taxon (1 .fasta for each locus)> -d <folder with contigs> -p <prefix in folder locus names relative to sequence names of added taxon> -n <label of added sequence> -t <number of threads>

## Define arguments 
while getopts r:d:p:n:a:s:t: opts
do
        case "${opts}"
		in
				r) added=${OPTARG};;
				d) dir=${OPTARG};;
				p) dirprefix=${OPTARG};;
				n) addname=${OPTARG};;
				a) alns=${OPTARG};;     # only needed if some reference sequences are in a separate alignment folder
				s) alnsuffix=${OPTARG};; # only needed if some reference sequences are in a separate alignment folder
				t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! $added ] ; then echo "folder with sequences to be added not specified (-r option), stopping" ; exit 0 ; fi
if [ ! $dir ] ; then echo "directory with exonerate results (-d option) not specified, stopping!" ; exit 0 ; fi
if [ ! $addname ] ; then echo "label of new sequence (-n option) not set, stopping" ; exit 0 ; fi

if [ ! -d $added ] ; then echo "folder with sequences to be added not found, stopping" ; exit 0 ; fi
if [ ! -d $dir ] ; then echo "directory with exonerate results <$dir> not found, stopping!" ; exit 0 ; fi

if [ ! $alns ] ; then echo "folder with alignment sequences to be added (backup sequences if absent in folder with sequences to be added) not specified (-a option), no backup search" ; alns='false' ; fi
if [ $alns != "false" ] && [ ! -d $alns ] ; then echo "folder with alignment sequences (backup sequences if absent in folder with sequences to be added) not found, stopping" ; exit 0 ; fi

if [ ! $alnsuffix ] ; then echo "suffix in alignment sequences not set (-s option) not set, setting to ''" ; alnsuffix='' ; fi
if [ ! $dirprefix ] ; then echo "prefix in exonerate locus name (-p option) not set, setting to ''" ; dirprefix='' ; fi

if [ ! $threads ] ; then echo "number of threads (-t option) not set, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
export dirsuffix=".bestScore.fasta"      # suffix in $dir folder
export addprefix=""                      # prefix in $added folder
export addsuffix=".fasta"                # suffix in $added folder
export alnprefix=""                      # prefix in $aln folder

export emptyseq="--------------------"   # string of sequences where no match could be found
export outprefix="add.contigs"           # prefix of outputted matching and non-matching loci

## Export variables
export added=$added
export dir=$dir
export dirprefix=$dirprefix
export dirsuffix=$dirsuffix
export addname=$addname
export alns=$alns
export alnprefix=$alnprefix
export alnsuffix=$alnsuffix

## Check that there are new sequences for each multifasta
echo "checking overlap between loci in ${dir} and ${added}"
if [ -f ${outprefix}.found ] ; then /bin/rm -f ${outprefix}.found ; fi
if [ -f ${outprefix}.lost ] ; then /bin/rm -f ${outprefix}.lost ; fi

example=$(ls -1d ${dir}/* | grep -v 'refseqs' | head -n 1)

for bfasta in $(ls -1d ${example}/*${dirsuffix})
do
	bfastabase=$(basename ${bfasta} ${dirsuffix})
	bsample=$(basename $(dirname $bfasta))

	addbase=$(echo ${bfastabase} | sed "s/^$bsample.$dirprefix//g")
	#addbase=$(echo ${bfastabase} | sed "s/^.*\($dirprefix.*\).*$/\1/" | sed "s/$dirprefix//g" # remove anything before $dirprefix, then remove leading $dirprefix
	toadd="${added}/${addprefix}${addbase}${addsuffix}"

	if [ -f $toadd ]
	then
		echo ${addbase} >> ${outprefix}.found
	else
		echo ${addbase} >> ${outprefix}.lost
	fi
done

# report results
if [ -f ${outprefix}.found ] ; then nfound=$(wc -l < ${outprefix}.found) ; else nfound=0 ; fi
if [ -f ${outprefix}.lost ] ; then nlost=$(wc -l < ${outprefix}.lost) ; else nlost=0 ; fi
nbfasta=$(ls -1d ${example}/${bsample}.${dirprefix}*${dirsuffix} 2> /dev/null | wc -l )
echo "found ${nbfasta} ${dirsuffix} files in ${example}"	
echo "found ${nfound} matching ${addsuffix} files in ${added}"


## Backup search
if [ $nlost -gt 0 ] && [ $alns != "false" ]
then
	echo "looking for ${nlost} sequences in ${alns}"
	nfound2=0

	for lost in $(cat ${outprefix}.lost)
	do
		aln=${alns}/${alnprefix}${lost}${alnsuffix}
		header=$(grep "^>${addname}$" ${aln})
		seq=$(grep "^>${addname}$" ${aln} -A 1 | grep -v '^>')
	
		if [ -z ${header} ]
		then
			echo ">${addname}" > ${added}/${addprefix}${lost}${addsuffix}
			echo ${emptyseq} >> ${added}/${addprefix}${lost}${addsuffix}
		else
			echo ${header} > ${added}/${addprefix}${lost}${addsuffix}
			echo ${seq} >> ${added}/${addprefix}${lost}${addsuffix}
			nfound2=$(expr ${nfound2} + 1)
		fi
	done
	
	# report results
	echo "found ${nfound2} matches in ${alns} (written to ${added})"	
fi

## Ask to continue
read -r -p "Continue? [y/n]" response
	case $response in
		[nN][oO]|[nN])
		exit 0
	esac

## Create empty sequences for each missing locus
if [ $nlost -gt 0 ] && [ $alns == "false" ]
then
	echo "creating empty files for ${nlost} sequences"
	for lost in $(cat ${outprefix}.lost)
	do
		echo ">${addname}" > ${added}/${addprefix}${lost}${addsuffix}
		echo ${emptyseq} >> ${added}/${addprefix}${lost}${addsuffix}
	done
fi

## Create new sample dir in $dir
if [ -d ${dir}/${addname} ] ; then echo "directory <${dir}/${addname}> already exists" ; exit 0 ; fi
mkdir ${dir}/${addname}

## Define function
doAdd()
{
	bfasta=$1
	
	bfastabase=$(basename ${bfasta} ${dirsuffix})
	bsample=$(basename $(dirname $bfasta))

	# get contig to be added
	addbase=$(echo ${bfastabase} | sed "s/^$bsample.$dirprefix//g")
	toadd="${added}/${addprefix}${addbase}${addsuffix}"
	addseq=$(cat ${toadd} | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | grep -v '^>')
	
	# add contig length to contig name
	clen=$(cat ${toadd} | awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' | grep -v '^>')
	cname="added_contig_${clen}_length"
	 
	fas="${dir}/${addname}/${addname}.${dirprefix}${addbase}${dirsuffix}"
	ascore="${dir}/${addname}/${addname}.${dirprefix}${addbase}.allScore"
	bscore="${dir}/${addname}/${addname}.${dirprefix}${addbase}.bestScore"
	efile="${dir}/${addname}/${addname}.${dirprefix}${addbase}.exonerate"
	
	echo $fas
	
	if [ -f $toadd ]
	then
		# .allScore
		echo "${cname} 5" > ${ascore} # assumes very low exonerate score!
		
		# .bestScore
		echo "${cname} 5" > ${bscore} # assumes very low exonerate score!
		
		# .exonerate
		touch $efile
		
		# .allScore.fasta
		echo ">${addname}" > ${fas}
		echo ${addseq} >> ${fas}
	else
		echo "sequence <$toadd> not found!" >> ${outprefix}.err
	fi

}

## Execute function
echo "adding sequences..."
export -f doAdd

ls -1d ${example}/${bsample}.${dirprefix}*${dirsuffix} | parallel -j $threads doAdd

## Finish
echo
echo "Done!"
echo
