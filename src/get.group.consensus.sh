#!/bin/bash

## Usage: get.group.consensus.sh -s <sample / group file> -d <indir> -m <minimum allele count / frequency> -b < minimum base count / frequency> -g <whether to ignore gaps> -n <whether to ignore ns> -a <prefix> -z <suffix> -v <visualize> -w <plot width> -h <plot height> -t <number of threads> 

## Needs: get.consensus.from.alignment.parallel.sh, create.multiconsensus.parallel.sh, align.multifastas.parallel.sh

## Define arguments
visualize='false'
ignoregaps='false'
ignoren='false'
while getopts s:d:m:b:gna:z:vw:h:t: opts
do
		case "${opts}"
		in
			s) sfile=${OPTARG};;
			d) indir=${OPTARG};;
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
if [ ! $sfile ] ; then echo "sample / group file not specified (-s option), stopping" ; exit 0 ; fi
if [ ! -f $sfile ] ; then echo "sample /group file <$sfile> not found, stopping" ; exit 0 ; fi

if [ ! $indir ] ; then echo "input directory with alignments not specified (-d option), stopping" ; exit 0 ; fi
if [ ! -d $indir ] ; then echo "input directory with alignments <$indir> not found, stopping" ; exit 0 ; fi

if [ ! ${minallfreq} ] ; then echo "minimum allele count / frequency (-m option) not specified, setting to 0.05", minallfreq=0.05 ; fi
if [ ! ${minbasefreq} ] ; then echo "minimum base count / frequency (-b option) not specified, setting to 0.5", minbasefreq=0.5 ; fi

if [ ${ignoregaps} == 'true' ] ; then echo "will ignore gaps to calculate minor allele frequencies (-g flag)" ; g="g" ; else g="" ; fi
if [ ${ignoren} == 'true' ] ; then echo "will omit N's from the consensus (-n flag)" ; n="n" ; else n="" ; fi
if [ ${visualize} == 'true' ] ; then echo "will visualize alignment and consensus (-v flag)" ; v="v" ; else v="" ; fi

if [ ! ${prefix} ] ; then echo "consensus sequence name prefix (-p option) not specified, setting to ''" ; prefix='' ; fi
if [ ! ${suffix} ] ; then echo "consensus sequence name suffix (-s option) not specified, setting to ''" ; suffix='' ; fi
if [ ! ${plotwidth} ] ; then echo "plot width (-w option) not specified, setting to 7" ; plotwidth=7 ; fi
if [ ! ${plotheight} ] ; then echo "plot height (-h option) not specified, setting to 15" ; plotheight=15 ; fi

if [ ! ${outdir} ] ; then echo "output directory (-o option) not specified, setting to '$(basename ${indir}).cons-${minallfreq}-${minbasefreq}'" ; outdir="$(basename ${indir}).cons-${minallfreq}-${minbasefreq}" ; fi

if [ ! ${threads} ] ; then echo "number of threads not specified (-t option), setting to 4" ; threads=4 ; fi
if [ ${threads} -gt 30 ] ; then echo "number of threads must not exceed 30, stopping" ; exit 0 ; fi

## Additional arguments
indirsuf=".fasta"
method="localpair"
argswitch="${g}${n}${v}"
argslen=$(expr length "$g$n$v")
if [ ${argslen} -gt 0 ] ; then args="-${argswitch}" ; else args="" ; fi

## Generate consensus per group
if [ -f grouppaths ] ; then rm grouppaths ; fi
for group in $(tail -n+2 ${sfile} | cut -f2 | sort | uniq | grep -v '^NA$')
do
	echo ${group}
	odir="${indir}.cons-${group}-${minallfreq}-${minbasefreq}"
	echo ${odir} >> grouppaths
	grep ${group} ${sfile} | cut -f1 > taxa_kept-${group}.txt
	get.consensus.from.alignment.parallel.sh -s taxa_kept-${group}.txt -d ${indir} -o ${odir} -m ${minallfreq} -b ${minbasefreq} -w ${plotwidth} -h ${plotheight} -t ${threads} ${args}
done

ngr=$(wc -l < grouppaths)
ls -1d $(head -n1 grouppaths)/*${indirsuf} | cut -f2 -d'/' | sed -e "s/${suffix}${indirsuf}//" > merged.loci.names
l=merged.loci.names
nloc=$(wc -l < $l)

## Align consensi per group
plen=$(expr length "${prefix}")
slen=$(expr length "${suffix}")
if [ ${plen} -gt 0 ] ; then pe="-p $prefix" ; else pe="" ; fi
if [ ${slen} -gt 0 ] ; then se="-s $suffix" ; else se="" ; fi

create.multiconsensus.parallel.sh -l ${l} -d grouppaths -t ${threads} ${pe} ${se} 
align.multifastas.parallel.sh -d multifasta.cons.${nloc}.${ngr} -m ${method} -t ${threads}

## Generate group consensus
get.consensus.from.alignment.parallel.sh -s grouppaths -d mafft.cons.${nloc}.${ngr} -m ${minallfreq} -b ${minbasefreq} -w ${plotwidth} -h ${plotheight} -t ${threads} ${args}

## Clean up
rm grouppaths
echo ; echo "Done" ; echo

