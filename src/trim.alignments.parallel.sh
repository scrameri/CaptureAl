#!/bin/bash

## Usage: trim.alignments.parallel.sh -s <sample / group file> -d <indir> -c <completeness for site trimming> -z <window size for sequence trimming> -n <fraction of required base matches in a window> -S <window step size> -i <turn on internal sequence trimming> -m <nachar> -v <turn on visualization> -p <turn on parsimony> -w <plotwidth> -h <plotheight> -t <threads>

## Needs: trim.alignment.R with ape library (and ggpubr library for parsimony plots)

## Define arguments
internal='false'
visualize='false'
parsimony='false'
while getopts s:d:c:z:n:S:im:vp:w:h:t: opts
do
        case "${opts}"
        in
        		s) sfile=${OPTARG};;
                d) indir=${OPTARG};;
                c) completeness=${OPTARG};;
                z) winsize=${OPTARG};;
                n) wincons=${OPTARG};;
                S) step=${OPTARG};;
                i) internal='false';; # FLAG
                m) nachar=${OPTARG};; 
                v) visualize='true';; # FLAG
                p) parsimony='true';; # FLAG
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

if [ ! $completeness ] ; then echo "completeness not specified (-c option), setting to 0.3" ; completeness=0.3 ; fi
if [ ! $winsize ] ; then echo "window size for slidign window approach not specified (-z option), setting to 20" ; winsize=20 ; fi
if [ ! $wincons ] ; then echo "fraction of base matches between window and consensus sequence not specified (-n option), setting to 0.5" ; wincons=0.5 ; fi
if [ ! $step ] ; then echo "sliding window step size not specified (-S option), setting to 1" ; step=1 ; fi
if [ $internal == 'true' ] ; then echo "will perform internal trimming (-i flag)" ; fi
if [ $visualize == 'true' ] ; then echo "will visualize trimming (-v flag)" ; fi
if [ $parsimony == 'true' ] ; then echo "will plot Parsimony Informative Sites (-p flag)" ; fi

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
export winsize=$winsize
export wincons=$wincons
export step=$step
export internal=$internal
export nachar=$nachar
export visualize=$visualize
export parsimony=$parsimony
export plotwidth=$plotwidth
export plotheight=$plotheight

## Create output directories
export odir="${indirbase}.c${completeness}.n${wincons}"
export logfolder="${odir}.logs"
if [ ! -d ${odir} ] ; then mkdir ${odir} ; else echo "${odir} already exists!" ; exit 0 ; fi
if [ ! -d ${logfolder} ] ; then mkdir ${logfolder} ; else echo "${logfolder} already exists!" ; exit 0 ; fi

## Create log file
echo "" >> ${logfolder}/params.log
echo "## LOGFILE for trim.alignments.parallel.sh ##" >> ${logfolder}/params.log
echo "" >> ${logfolder}/params.log
echo "sample / group file:     ${sfile}" >> ${logfolder}/params.log
echo "input directory:         ${indir}" >> ${logfolder}/params.log
echo "completeness:            ${completeness}" >> ${logfolder}/params.log
echo "winsize:                 ${winsize}" >> ${logfolder}/params.log
echo "wincons:                 ${wincons}" >> ${logfolder}/params.log
echo "step:                    ${step}" >> ${logfolder}/params.log
echo "internal:                ${internal}" >> ${logfolder}/params.log
echo "NA.char (encoding gaps): ${nachar}" >> ${logfolder}/params.log
echo "visualize:               ${visualize}" >> ${logfolder}/params.log
echo "parsimony:               ${parsimony}" >> ${logfolder}/params.log
echo "plot width:              ${plotwidth}" >> ${logfolder}/params.log
echo "plot height:             ${plotheight}" >> ${logfolder}/params.log
echo "threads:                 ${threads}" >> ${logfolder}/params.log

## Define function
doTrim()
{
	aln=$1
	alnbase=$(basename $aln ${alnsuffix})
	echo $alnbase

	args="$sfile $aln $odir $completeness $winsize $wincons $step $internal $nachar $visualize $parsimony $plotwidth $plotheight"
	#echo $args

	trim.alignment.R ${args} 1> ${logfolder}/${alnbase}.itr.log 2> ${logfolder}/${alnbase}.itr.err

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
