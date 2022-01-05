#!/bin/bash

## Usage: add.seq.to.multifasta.sh -r <folder with sequences of 1 taxon (1 .fasta for each locus)> -m <folder with multifastas> -p <prefix in multifasta sequence names relative to sequence names of added taxon> -n <label of added sequence> -a <dir with aligned sequences for backup search> -s <suffix of aligned sequences> -t <number of threads>

## Define arguments 
while getopts r:m:p:n:a:s:t: opts
do
        case "${opts}"
		in
				r) added=${OPTARG};;
				m) mfastas=${OPTARG};;
				p) mfastaprefix=${OPTARG};;
				n) addname=${OPTARG};;
				a) alns=${OPTARG};;     # only needed if some reference sequences are in a separate alignment folder
				s) alnsuffix=${OPTARG};; # only needed if some reference sequences are in a separate alignment folder
				t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! $added ] ; then echo "folder with sequences to be added not specified (-r option), stopping" ; exit 0 ; fi
if [ ! $mfastas ] ; then echo "folder with multifasta files (where sequences should be added to) not specified (-m option), stopping" ; exit 0 ; fi
if [ ! $alns ] ; then echo "folder with alignment sequences to be added (backup sequences if absent in folder with sequences to be added) not specified (-a option), no backup search" ; alns='false' ; fi

if [ ! -d $added ] ; then echo "folder with sequences to be added not found, stopping" ; exit 0 ; fi
if [ ! -d $mfastas ] ; then echo "folder with multifasta files (where sequences should be added to) not found, stopping" ; exit 0 ; fi
if [ $alns != "false" ] && [ ! -d $alns ] ; then echo "folder with alignment sequences (backup sequences if absent in folder with sequences to be added) not found, stopping" ; exit 0 ; fi

if [ ! $mfastaprefix ] ; then echo "prefix in multifasta sequence name (-p option) not set, setting to ''" ; mfastaprefix='' ; fi
if [ ! $addname ] ; then echo "label of new sequence (-n option) not set, stopping" ; exit 0 ; fi
if [ ! $alnsuffix ] ; then echo "suffix in alignment sequences not set (-s option) not set, setting to ''" ; alnsuffix='' ; fi

if [ ! $threads ] ; then echo "number of threads (-t option) not set, setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
export addsuffix=".fasta"                # suffix in $added folder
export mfastasuffix=".all.fasta"         # suffix in $mfastas folder
export emptyseq="--------------------"   # string of sequences where no match could be found
export outprefix="add.seq.to.multifasta" # prefix of outputted matching and non-matching loci

## Export variables
export added=$added
export alns=$alns
export alnsuffix=$alnsuffix
export mfastas=$mfastas
export mfastaprefix=$mfastaprefix
export addname=$addname

## Check that there are new sequences for each multifasta
echo "checking overlap between loci in ${mfastas} and ${added}"
if [ -f ${outprefix}.found ] ; then /bin/rm -f ${outprefix}.found ; fi
if [ -f ${outprefix}.lost ] ; then /bin/rm -f ${outprefix}.lost ; fi

for mfasta in $(ls -1d ${mfastas}/*${mfastasuffix})
	do
		mfastabase=$(basename ${mfasta} ${mfastasuffix})
		addbase=$(echo ${mfastabase} | sed -e s/${mfastaprefix}//g)
		toadd="${added}/${addbase}${addsuffix}"
	
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
naln=$(ls -1d ${mfastas}/*${mfastasuffix} | wc -l)
echo "found ${naln} multifasta files in ${mfastas}"	
echo "found ${nfound} matching fasta files in ${added}"


## Concatenate ${outprefix}.lost with a spacer
if [ $nlost -gt 0 ] && [ $alns != "false" ]
then
	echo "looking for ${nlost} sequences in ${alns}"
	nfound2=0

	for lost in $(cat ${outprefix}.lost)
	do
		aln=${alns}/${lost}${alnsuffix}
		header=$(grep "^>${addname}$" ${aln})
		seq=$(grep "^>${addname}$" ${aln} -A 1 | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | grep -v '^>')
	
		if [ -z ${header} ]
		then
			echo ">${addname}" > ${added}/${lost}${addsuffix}
			echo ${emptyseq} >> ${added}/${lost}${addsuffix}
		else
			echo ${header} > ${added}/${lost}${addsuffix}
			echo ${seq} >> ${added}/${lost}${addsuffix}
			nfound2=$(expr ${nfound2} + 1)
		fi
	done
	
	# report results
	echo "found ${nfound2} matches in ${alns} (written to ${added})"
	
fi
if [ $nlost -gt 0 ] && [ $alns == "false" ]
then
	echo "creating empty files for ${nlost} sequences"
	for lost in $(cat ${outprefix}.lost)
	do
		echo ">${addname}" > ${added}/${lost}${addsuffix}
		echo ${emptyseq} >> ${added}/${lost}${addsuffix}
	done
fi

## Ask to continue
read -r -p "Continue? [y/n]" response
	case $response in
		[nN][oO]|[nN])
		exit 0
	esac

## Define function
doAdd()
{
	mfasta=$1
	mfastabase=$(basename $mfasta $mfastasuffix)
	echo ${mfastabase}
	
	addbase=$(echo ${mfastabase} | sed -e s/${mfastaprefix}//g)
	toadd="${added}/${addbase}${addsuffix}"
	addseq=$(cat ${toadd} | awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' | grep -v '^>')
	
	if [ -f $toadd ]
	then
		ispresent=$(grep "^>${addname}$" ${mfasta})
		if [ -z ${ispresent} ]
		then
			echo ">${addname}" >> ${mfasta}
			echo ${addseq} >> ${mfasta}
		else 
			echo "sequence <$addname> already in ${mfasta}, not added again"
		fi
	fi
}

## Execute function
echo "appending new sequences..."
export -f doAdd

ls -1d ${mfastas}/*${mfastasuffix} | parallel -j $threads doAdd

## Finish
echo
echo "done!"
echo
