#!/bin/bash

## Usage: filter.visual.assemblies.sh -s <taxon group file> -t <loci_stats.txt> -r <refseqs.fasta> -a <minpregion> -b <minpcontig> -c <minptaxa> -d <maxncontigs> -e <minbestscorenorm> -f <mintaln> -g <mintfrac> -h <minbestscore> -i <minbestlength> -p <minfrac>

## Description
# Use up to 10 filtering criteria to identify target regions with adequate assembly quality across all defined taxon groups. 

## Details
# The three required arguments are:
# -s	<sfile>			the file mapping samples to taxon groups
# -t	<stats>			the EXONERATE statistics
# -r	<refseqs>		the target region reference sequences (FASTA)
# 
# The first three filters take absolute thresholds and aim to remove poorly assembled samples or target regions:
# -a	<minpregion>	minimum fraction of regions with at least one contig in a sample (filters samples)
# -b	<minpcontig>	minimum median contig length relative to target length (filters samples)
# -c	<minptaxa>		minimum fraction of samples with at least one contig in a region (filters target regions)
# 
# The next six filters take thresholds that need to be met in a specified fraction of samples in each considered taxon group:
# -d	<maxncontigs>	maximum number of non-zero (fragments combined) contigs in a target region
# -e	<minbestscorenorm>	minimum normalized EXONERATE alignment score
# -f	<mintaln>		minimum EXONERATE alignment length
# -g	<mintfrac>		minimum alignment fraction (EXONERATE alignment length divided by target region length)
# -h	<minbestscore>	minimum raw EXONERATE alignment score
# -i	<minbestlength>	minimum length of best-matching contig
# 
# The last filter regulates how strictly the six filters above are applied
# -p	<minfrac>		minimum fraction of samples in each taxon group that need to pass each filter in order to keep a certain target region


## Define arguments
while getopts s:t:r:a:b:c:d:e:f:g:h:i:p: opts
do
		case "${opts}"
		in
			s) sfile=${OPTARG};;
			t) stats=${OPTARG};;
			r) refseqs=${OPTARG};;
			a) minpregion=${OPTARG};;
			b) minpcontig=${OPTARG};;
			c) minptaxa=${OPTARG};;
			d) maxncontigs=${OPTARG};;
			e) minbestscorenorm=${OPTARG};;
			f) mintaln=${OPTARG};;
			g) mintfrac=${OPTARG};;
			h) minbestscore=${OPTARG};;
			i) minbestlength=${OPTARG};;
			p) minfrac=${OPTARG};;
		esac
done

## Default parameters
# required
if [ ! $sfile ] ; then echo "<sfile> parameter not specified (-s option), stopping" ; exit 0 ; fi
if [ ! $stats ] ; then echo "<stats> parameter not specified (-t option), stopping" ; exit 0 ; fi
if [ ! $refseqs ] ; then echo "<refseqs> parameter not specified (-r option), stopping" ; exit 0 ; fi

# optional
if [ ! $minpregion ] ; then echo "<minpregion> filter not set (-a option), setting minpregion=0.3" ; minpregion=0.3 ; fi
if [ ! $minpcontig ] ; then echo "<minpcontig> filter not set (-b option), setting minpcontig=0.5" ; minpcontig=0.5 ; fi
if [ ! $minptaxa ] ; then echo "<minptaxa> filter not set (-c option), setting minptaxa=0.3" ; minptaxa=0.3 ; fi

if [ ! $maxncontigs ] ; then echo "<maxncontigs> filter not set (-d option), setting maxncontigs=1" ; maxncontigs=1 ; fi
if [ ! $minbestscorenorm ] ; then echo "<minbestscorenorm> filter not set (-e option), setting minbestscorenorm=1" ; minbestscorenorm=1 ; fi
if [ ! $mintaln ] ; then echo "<mintaln> filter not set (-f option), setting mintaln=50" ; mintaln=50 ; fi
if [ ! $mintfrac ] ; then echo "<mintfrac> filter not set (-g option), setting mintfrac=0.2" ; mintfrac=0.2 ; fi
if [ ! $minbestscore ] ; then echo "<minbestscore> filter not set (-h option), setting minbestscore=1" ; minbestscore=1 ; fi
if [ ! $minbestlength ] ; then echo "<minbestlength> filter not set (-i option), setting minbestlength=1" ; minbestlength=1 ; fi

if [ ! $minfrac ] ; then echo "<minfrac> filter not set (-p option), setting minfrac=0.9" ; minfrac=0.9 ; fi

## Execute R script
filter.visual.assemblies.R ${sfile} ${stats} ${ref} ${minpregion} ${minpcontig} ${minptaxa} ${maxncontigs} ${minbestscorenorm} ${mintaln} ${mintfrac} ${minbestscore} ${minbestlength} ${minfrac}
