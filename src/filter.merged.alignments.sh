#!/bin/bash

## Usage: filter.merged.alignments.sh -d <directory with merged alignments> -s <score>

## Author: simon.crameri@env.ethz.ch, Apr 2019

## Define arguments
while getopts d:s: opts
do
        case "${opts}"
        in
                d) dir=${OPTARG};;  
                s) mthr=${OPTARG};;
         esac
done

## Check arguments
if [ ! $dir ] ; then echo "directory with merged alignments not specified (-d option)" ; exit 0 ; fi
if [ ! -d $dir ] ; then echo "directory with merged alignments not found" ; exit 0 ; fi

if [ ! $mthr ] ; then echo "alignment score threshold not specified (-s option), setting to 0.95" ; mthr=0.95 ; fi
mthrtest=$(echo $mthr | awk '{ if ($1<0 || $1>1) { print "false" } else { print "true" } }')
if [ $mthrtest == 'false' ] ; then echo "alignment score ($mthr) must be between 0 and 1" ; exit 0 ; fi

## Additional arguments
odir="${dir}/loci_rm_score_${mthr}"

## Get alignments to be filtered
if [ -f ${dir}.score_${mthr} ] ; then mv ${dir}.score_${mthr} ${dir}.score_${mthr}.bak ; fi
cat ${dir}/*log | grep 'Overall alignment score:' | awk -F '\t' -v thr=$mthr '($2 < thr)' | cut -f3 | cut -f2 -d'/' > ${dir}.score_${mthr}
nfilt=$(wc -l < ${dir}.score_${mthr})

## Reverse any previous filterings
nrmfolder=$(ls -1d ${dir}/loci_rm_score_* 2> /dev/null | wc -l)
if [ $nrmfolder -gt 0 ]
then
	for rmfolder in $(ls -1d ${dir}/loci_rm_score_*)
	do
		mv ${rmfolder}/* ${dir}/ 2> /dev/null
		/bin/rm -r ${rmfolder}
	done
fi

## Create output dir
if [ $nfilt -gt 0 ] ; then mkdir ${odir} ; fi

## Move filtered alignments
for path in $(cat ${dir}.score_${mthr}) 
do
	filename=$(basename "$path")
	ext=".${filename##*.}"
	pathbase="${filename%.*}"
	
	mv ${dir}/${pathbase}.fasta ${odir}
	# mv ${dir}/${pathbase}.log ${odir} # log files stay in $dir!
	mv ${dir}/${pathbase}.pdf ${odir} 
done

## Finish
if [ $nfilt -eq 0 ] ; then /bin/rm -f ${dir}.score_${mthr} ; fi
echo 
echo "filtered ${nfilt} alignment(s) with score below ${mthr}!"
echo