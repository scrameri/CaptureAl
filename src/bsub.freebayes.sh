#!/bin/bash

## Resource usage
#BSUB -J "SNP[1-277]%12"              ## freebayes runs [1-X] in %Y parallel jobs
#BSUB -W 12:00                            ## maximum runtime of X:00 hours for each job 
#BSUB -n 1                                ## n core(s) per job -> n * Y cores in parallel
#BSUB -R "rusage[mem=10000]"              ## Faster queue if total R (R*%Y) < 130000 (100000 in reduced mode over holidays)

## Memory notes
# before running freebayes, you might want to cap .bam files to a maximum coverage of e.g. 500 (run bsub.bwamem.sh with maxcov=500 or as low as 50)
# if memory is still an issue, set --use-best-n-alleles to 4 to 7  (see https://github.com/freebayes/freebayes/issues/582)
# if memory is still an issue, use a merged .bam file (run submit.merge.bamfiles.sh and pass merged.bam to freebayes via -b option) (N. Zemp, 29 Dec. 2021)
# if memory is still an issue, use --limit-coverage to 20 or 30 maximum coverage (see https://github.com/freebayes/freebayes/issues/582)
# if memory is still an issue, set --skip-coverage to some acceptable number, e.g. 100 x the number of samples (see https://github.com/freebayes/freebayes/issues/582)

# arguments
lenperjob=10000
#lenperjob=7500
bedfolder="regions_${lenperjob}"
ofolder="vcfs"
#bamlist="bamlist.txt" # for -L option (uses more memory)
mergedbam="merged.bam" # for -b option (needs merged .bam file, much more memory-efficient)
ref="consDalbergia_4c_2396.fasta"
#nsamples=$(cat ${bamlist} |wc -l)
#maxcov=$((500*${nsamples})) # for --skip-coverage
maxcov=20 # for --limit-coverage

# check arguments
if [ ! -d ${ofolder} ] ; then mkdir ${ofolder} ; fi
if [ ! -d ${bedfolder} ] ; then echo "input directory (.bam files) <${bedfolder}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${ref} ] ; then echo "reference .fasta file <${ref}> not found, stopping" ; exit 0 ; fi
if [ ! -f ${mergedbam} ] ; then echo "merged .bam file <${mergedbam}> not found, stopping" ; exit 0 ; fi
#if [ ! -f ${bamlist} ] ; then echo "list of .bam files <${bamlist}> not found, stopping" ; exit 0 ; fi

# job (region) index
IDX=$LSB_JOBINDEX

# load modules
module load gcc/4.8.2 gdc python/2.7.11 freebayes/1.3.4

# run freebayes
# options in v. 1.3.4:
#-L --bam-list FILE				A file containing a list of BAM files to be analyzed.
#-b --bam FILE  					Add FILE to the set of BAM files to be analyzed.
#-f --fasta-reference 			FILE	Use FILE as the reference sequence for analysis. An index file (FILE.fai) will be created if none exists. If neither --targets nor --region are specified, FreeBayes will analyze every position in this reference.
#-t --targets FILE				Limit analysis to targets listed in the BED-format FILE.
#-p --ploidy N   				Sets the default ploidy for the analysis to N.  default: 2
#-C --min-alternate-count N		Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position.  default: 2
#-F --min-alternate-fraction N	Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position.  default: 0.05
#--limit-coverage N				Downsample per-sample coverage to this level if greater than this coverage. default: no limit
#-g --skip-coverage N			Skip processing of alignments overlapping positions with coverage >N. This filters sites above this coverage, but will also reduce data nearby.  default: no limit # Github explanation: This is for the total coverage among all samples. If coverage goes above this, the reads are only counted until coverage goes back below the threshold.
#-n --use-best-n-alleles N		Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores.  (Set to 0 to use all; default: all)
#-E --max-complex-gap N --haplotype-length N	Allow haplotype calls with contiguous embedded matches of up to this length. Set N=-1 to disable clumping. (default: 3)
#-v --vcf FILE   				Output VCF-format results to FILE. (default: stdout)

# for freebayes/1.1.0-3-g961e5f3
#freebayes -f $ref -L $bamlist -t ${bedfolder}/regions_${lenperjob}_${IDX} -v ${ofolder}/raw.${IDX}.vcf --min-alternate-fraction 0.05 --min-repeat-entropy 1 --use-best-n-alleles 4

# for freebayes/1.3.4 using -L bamlist.txt
#freebayes -L ${bamlist} -f $ref -t ${bedfolder}/regions_${lenperjob}_${IDX} \
#          -p 2 -C 1 -F 0.05 --skip-coverage ${maxcov} \
#          --use-best-n-alleles 4 --max-complex-gap -1 --haplotype-length -1  \
#          -v ${ofolder}/raw.${IDX}.vcf

# for freebayes/1.3.4 using -b merged.bam file
freebayes -b ${mergedbam} -f $ref -t ${bedfolder}/regions_${lenperjob}_${IDX} \
          -p 2 -C 1 -F 0.05 \
          --use-best-n-alleles 4 -E -1 \
          -v ${ofolder}/raw.${IDX}.vcf

