#!/bin/bash

## Usage: filter.visual.coverages.sh -s <taxon group file> -t <coverage_stats.txt> -r <refseqs.fasta> -a <minpregion> -b <minptaxa> -c <minlen> -d <mincov> -e <maxcov> -f <minratio> -p <minfrac>

## Description
# Use up to 7 filtering criteria to identify target regions with adequate sequencing data across all defined taxon groups. 

## Details
# The three required arguments are:
# -s	<sfile>			the file mapping samples to taxon groups
# -t	<stats>			the COVERAGE statistics
# -r	<refseqs>		the target region reference sequences (FASTA)
# 
# The first two filters take absolute thresholds and aim to remove poorly sequenced samples or target regions:
# -a	<minpregion>	minimum fraction of regions with at least one mapped read in a sample (filters samples)
# -b	<minptaxa>		minimum fraction of samples with at least one mapped read in a region (filters target regions)
# 
# The next four filters take thresholds that need to be met in a specified fraction of samples in each considered taxon group:
# -c	<minlen>		minimum BWA-MEM alignment length
# -d	<mincov>		minimum average coverage in the aligned region
# -e	<maxcov>		maximum average coverage in the aligned region
# -f	<minratio>		minimum alignment fraction (BWA-MEM alignment length divided by target region length)
# 
# The last filter regulates how strictly the four filters above are applied
# -p	<minfrac>		minimum fraction of samples in each taxon group that need to pass each filter in order to keep a certain target region


## Define arguments
while getopts s:t:r:a:b:c:d:e:f:p: opts
do
		case "${opts}"
		in
			s) sfile=${OPTARG};;
			t) stats=${OPTARG};;
			r) refseqs=${OPTARG};;
			a) minpregion=${OPTARG};;
			b) minptaxa=${OPTARG};;
			c) minlen=${OPTARG};;
			d) mincov=${OPTARG};;
			e) maxcov=${OPTARG};;
			f) minratio=${OPTARG};;
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
if [ ! $minptaxa ] ; then echo "<minptaxa> filter not set (-b option), setting minptaxa=0.3" ; minptaxa=0.3 ; fi

if [ ! $minlen ] ; then echo "<minlen> filter not set (-c option), setting minlen=500" ; minlen=500 ; fi
if [ ! $mincov ] ; then echo "<mincov> filter not set (-d option), setting mincov=10" ; mincov=10 ; fi
if [ ! $maxcov ] ; then echo "<maxcov> filter not set (-e option), setting maxcov=1000" ; maxcov=1000 ; fi
if [ ! $minratio ] ; then echo "<minratio> filter not set (-f option), setting minratio=0.5" ; minratio=0.5 ; fi

if [ ! $minfrac ] ; then echo "<minfrac> filter not set (-p option), setting minfrac=0.9" ; minfrac=0.9 ; fi

## Execute R script
filter.visual.coverages.R ${sfile} ${stats} ${ref} ${minpregion} ${minptaxa} ${minlen} ${mincov} ${maxcov} ${minratio} ${minfrac} 
