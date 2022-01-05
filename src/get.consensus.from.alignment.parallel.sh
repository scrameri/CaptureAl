#!/bin/bash

## Usage: get.consensus.from.alignment.sh -s <sample file> -d <dir with alignments> -m <minimum allele count / frequency> -b < minimum base count / frequency> -g <whether to ignore gaps> -n <whether to ignore ns> -a <prefix> -z <suffix> -v <visualize> -w <plot width> -h <plot height> -t <number of threads> 

## Needs:
# get.consensus.from.alignment.R (needs R ape library)

## Define arguments
visualize='false'
ignoregaps='false'
ignoren='false'
while getopts s:d:o:m:b:gna:z:vw:h:t: opts
do
	case "${opts}"
		in
			s) sfile=${OPTARG};;
			d) indir=${OPTARG};;
			o) outdir=${OPTARG};;
			m) minallfreq=${OPTARG};;
			b) minbasefreq=${OPTARG};;
			g) ignoregaps='true';;
			n) ignoren='true';;
			a) prefix=${OPTARG};;    # use 'FALSE' instead of an empty string
			z) suffix=${OPTARG};;    # use 'FALSE' instead of an empty string
			v) visualize='true';; 
			w) plotwidth=${OPTARG};;
			h) plotheight=${OPTARG};;
			t) threads=${OPTARG};;
		esac
done

## Check arguments
if [ ! ${sfile} ] ; then echo "sample file (-s option) not specified" ; exit 0 ; fi
if [ ! -f ${sfile} ] ; then echo "sample file  <${sfile}> does not exist" ; exit 0 ; fi

if [ ! ${indir} ] ; then echo "input directory (-d option) not specified" ; exit 0 ; fi
if [ ! -d ${indir} ] ; then echo "input directory <${indir}> does not exist" ; exit 0 ; fi

if [ ! ${minallfreq} ] ; then echo "minimum allele count / frequency (-m option) not specified, setting to 0.05", minallfreq=0.05 ; fi
if [ ! ${minbasefreq} ] ; then echo "minimum base count / frequency (-b option) not specified, setting to 0.5", minbasefreq=0.5 ; fi

if [ ${ignoregaps} == 'true' ] ; then echo "will ignore gaps to calculate minor allele frequencies (-g flag)" ; fi
if [ ${ignoren} == 'true' ] ; then echo "will omit N's from the consensus (-n flag)" ; fi
if [ ${visualize} == 'true' ] ; then echo "will visualize alignment and consensus (-v flag)" ; fi

if [ ! ${prefix} ] ; then echo "consensus sequence name prefix (-p option) not specified, setting to ''" ; prefix='FALSE' ; fi
if [ ! ${suffix} ] ; then echo "consensus sequence name suffix (-s option) not specified, setting to ''" ; suffix='FALSE' ; fi
if [ ! ${plotwidth} ] ; then echo "plot width (-w option) not specified, setting to 7" ; plotwidth=7 ; fi
if [ ! ${plotheight} ] ; then echo "plot height (-h option) not specified, setting to 15" ; plotheight=15 ; fi

if [ ! ${outdir} ] ; then echo "output directory (-o option) not specified, setting to '$(basename ${indir}).cons-${minallfreq}-${minbasefreq}'" ; outdir="$(basename ${indir}).cons-${minallfreq}-${minbasefreq}" ; fi

if [ ! ${threads} ] ; then echo "number of threads (-t option) not specified, setting to 4" ; threads=4 ; fi
if [ ${threads} -gt 30 ] ; then echo "number of threads must not exceed 30" ; exit 0 ; fi

## Additional arguments
export alnsuffix=".fasta"
#outdir="$(basename ${indir}).cons-${minallfreq}-${minbasefreq}"

## Create output dirs and move any existing dirs or files
if [ -d ${outdir} ] ; then mv ${outdir} ${outdir}.bak ; fi
mkdir ${outdir}

logfolder=${outdir}.logs
if [ -d ${logfolder} ] ; then mv ${logfolder} ${outdir}.bak.logs ; fi
mkdir ${logfolder}

vizfolder=${outdir}.viz
if [ -d ${vizfolder} ] ; then mv ${vizfolder} ${outdir}.bak.viz ; fi

if [ -f ${outdir}${alnsuffix} ] ; then mv ${outdir}${alnsuffix} ${outdir}.bak ; fi

## Export arguments
export sfile=$sfile
export indir=$indir
export outdir=$outdir
export outdir=$outdir
export logfolder=$logfolder
export minallfreq=$minallfreq
export minbasefreq=$minbasefreq
export ignoregaps=$ignoregaps
export ignoren=$ignoren
export prefix=$prefix
export suffix=$suffix
export visualize=$visualize
export plotwidth=$plotwidth
export plotheight=$plotheight

## Create log file
echo "" > ${logfolder}/params.log
echo "## LOGFILE for get.consensus.from.alignment.parallel.sh ##" >> ${logfolder}/params.log
echo "" >> ${logfolder}/params.log
echo "sample file:                            ${sfile}" >> ${logfolder}/params.log
echo "input directory:                        ${indir}" >> ${logfolder}/params.log
echo "output directory:                       ${outdir}" >> ${logfolder}/params.log
echo "minimum allele freq/count:              ${minallfreq}" >> ${logfolder}/params.log
echo "minimum base freq/count:                ${minbasefreq}" >> ${logfolder}/params.log
echo "ignore gaps when computing allele freq: ${ignoregaps}" >> ${logfolder}/params.log
echo "omit N from consensus:                  ${ignoren}" >> ${logfolder}/params.log
echo "prefix:                                 ${prefix}" >> ${logfolder}/params.log
echo "suffix:                                 ${suffix}" >> ${logfolder}/params.log
echo "visualize:                              ${visualize}" >> ${logfolder}/params.log
echo "plot width:                             ${plotwidth}" >> ${logfolder}/params.log
echo "plot height:                            ${plotheight}" >> ${logfolder}/params.log
echo "threads:                                ${threads}" >> ${logfolder}/params.log

## Define function
make_cons()
{
	aln=$1
	alnbase=$(basename $aln ${alnsuffix})
	echo $alnbase

	args="$aln $sfile $outdir $minallfreq $minbasefreq $ignoregaps $ignoren $prefix $suffix $visualize $plotwidth $plotheight"
	get.consensus.from.alignment.R $args 1> ${logfolder}/${alnbase}.cons.log 2> ${logfolder}/${alnbase}.cons.err
}

## Execute function
export -f make_cons
ls -1d ${indir}/*${alnsuffix} | parallel -j ${threads} make_cons

## Concatenate .fastas
echo 
echo "Concatenating consensus ${alnsuffix} files..."
cat ${outdir}/*.cons${alnsuffix} > ${outdir}${alnsuffix}

## Finish
echo 
echo "All sequences processed!"
echo
