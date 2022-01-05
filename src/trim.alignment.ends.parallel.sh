#!/bin/bash

## Usage: trim.alignment.ends.parallel.sh -s <sample / group file> -d <indir> -c <completeness> -n <maxnucdiv> -m <nachar> -v <flag for visualization> -w <plotwidth> -h <plotheight> -t <threads>

## Needs: trim.alignment.ends.R with ape library

## Define arguments
visualize='false'
while getopts s:d:c:n:m:vw:h:t: opts
do
        case "${opts}"
        in
        		s) sfile=${OPTARG};;
                d) indir=${OPTARG};;
                c) completeness=${OPTARG};;
                n) maxnucdiv=${OPTARG};;
                m) nachar=${OPTARG};;
                v) visualize='true';;
                w) plotwidth=${OPTARG};;
                h) plotheight=${OPTARG};;
                t) threads=${OPTARG};;
         esac
done

## Check arguments
if [ ! $sfile ] ; then echo "sample / group file not specified (-s option), stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample /group file <$sfile> not found, stopping" ; exit 0 ; fi

if [ ! $indir ] ; then echo "input directory with alignments not specified (-d option), stopping" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "input directory with alignments <$indir> not found, stopping" ; exit 0 ; fi

if [ ! $completeness ] ; then echo "completeness not specified (-c option), setting to 0.5" ; completeness=0.5 ; fi
if [ ! $maxnucdiv ] ; then echo "maximum nucleotide diversity not specified (-n option), setting to 0.25" ; maxnucdiv=0.25 ; fi
if [ $visualize == 'true' ] ; then echo "will visualize trimming (-v flag)" ; fi

if [ ! $nachar ] ; then nachar='-' ; fi
if [ ! $plotwidth ] ; then plotwidth=15 ; fi
if [ ! $plotheight ] ; then plotheight=7 ; fi

if [ ! $threads ] ; then echo "number of threads not specified (-t option), setting to 4" ; threads=4 ; fi
if [ $threads -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
export alnsuffix=".fasta"
export indirbase=$(basename $indir)

## Export arguments
export sfile=$sfile
export indir=$indir
export completeness=$completeness
export maxnucdiv=$maxnucdiv
export nachar=$nachar
export visualize=$visualize
export plotwidth=$plotwidth
export plotheight=$plotheight

## Create output directories
export odir="${indirbase}.c${completeness}.d${maxnucdiv}"
export logfolder="${odir}.logs"
if [ ! -d ${odir} ] ; then mkdir ${odir} ; else echo "${odir} already exists!" ; exit 0 ; fi
if [ ! -d ${logfolder} ] ; then mkdir ${logfolder} ; else echo "${logfolder} already exists!" ; exit 0 ; fi

## Create log file
echo "" >> ${logfolder}/params.log
echo "## LOGFILE for trim.alignment.ends.parallel.sh ##" >> ${logfolder}/params.log
echo "" >> ${logfolder}/params.log
echo "sample / group file:       ${sfile}" >> ${logfolder}/params.log
echo "input directory:           ${indir}" >> ${logfolder}/params.log
echo "completeness:              ${completeness}" >> ${logfolder}/params.log
echo "max. nucleotide diversity: ${maxnucdiv}" >> ${logfolder}/params.log
echo "NA.char (encoding gaps):   ${nachar}" >> ${logfolder}/params.log
echo "visualize:                 ${visualize}" >> ${logfolder}/params.log
echo "plot width:                ${plotwidth}" >> ${logfolder}/params.log
echo "plot height:               ${plotheight}" >> ${logfolder}/params.log
echo "threads:                   ${threads}" >> ${logfolder}/params.log

## Define function
doTrim()
{
	aln=$1
	alnbase=$(basename $aln ${alnsuffix})
	echo $alnbase

	args="$sfile $aln $odir $completeness $maxnucdiv $nachar $visualize $plotwidth $plotheight"
	#echo $args

	trim.alignment.ends.R ${args} 1> ${logfolder}/${alnbase}.etr.log 2> ${logfolder}/${alnbase}.etr.err

}
export -f doTrim

## Execute function
ls -1 ${indir}/*${alnsuffix} | parallel -j ${threads} doTrim

## Be verbose
grep -l 'kept 0 /' ${logfolder}/*log | cut -f2 -d'/' | sed s/.log$/${alnsuffix}/g > ${odir}.zeroloci
nfiltered=$(wc -l < ${odir}.zeroloci)
echo ; echo "$nfiltered loci were completely trimmed!"
if [ $nfiltered -eq 0 ] ; then /bin/rm -f ${odir}.zeroloci ; fi

## Finish
echo
echo "All alignments processed!"
echo
