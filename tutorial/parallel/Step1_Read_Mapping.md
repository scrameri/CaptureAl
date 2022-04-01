[← Back to Read Trimming](Step0.2_Read_Trimming.md)


# STEP 1

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step1.png)


## 1) [run.bwamem.sh](https://github.com/scrameri/CaptureAl/wiki/run.bwamem.sh)

This maps the reads of samples specified in [samples.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/samples.txt) located in the `NovaSeq_run1_trimmed` directory against a [reference.fasta](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/reference.fasta) FASTA file.

It assumes *paired-end* reads located in separate files, e.g. `BEN001.trim1.fastq.gz` and `BEN001.trim2.fastq.gz` for sample *BEN001*. It further processes 4 samples in parallel, but performs hyper-threading for each sample, using 5 times as many threads as specified in `-t`.

After mapping using a `-T` quality threshold, the BAM files are further filtered for high-quality mappings (`-Q` option), and PCR duplicates are removed using [Picard Tools](https://broadinstitute.github.io/picard/).

```
run.bwamem.sh -s samples.txt -r reference.fasta -e '.trim1.fastq.gz,.trim2.fastq.gz' -T 10 -Q 20 \
              -d NovaSeq_run1_trimmed -t 4
```

## 2) [get.coverage.stats.sh](https://github.com/scrameri/CaptureAl/wiki/get.coverage.stats.sh)

**Usage**
```
get.coverage.stats.sh -Q <positive integer> -d <directory> -s <file> -t <integer>
```

**Arguments**
```
# Required
-Q                minimum mapping quality used when running run.bwamem.sh 


# Optional [DEFAULT]
-d  [pwd]         path to directory with mapping results (directory with sample subdirectories containing mapped reads)
-s  [samples.txt] file with sample names (without header or '>')
-t  [2]           number of samples to process in parallel

```

**Depends on**
```
get.coverage.stats.R
samtools
bedtools
```


**Example**
```
get.coverage.stats.sh -s samples.txt -Q 20 -t 20
```

## 3) collect.coverage.stats.R

**Usage**
```
collect.coverage.stats.R <file> <positive integer>
```

**Arguments**
```
# Required
1) <sfile|CHR>:  path to sample file containing paths to sample directories with mapping results. No header expected. ;

# Optional [DEFAULT] (if one or more arguments are passed, they must be passed in this order)
2) <Q|NUM>:      mapping quality threshold. Will look for SAMPLE.Q${Q}.coverage.txt [DEFAULT: 10]
3) <maxcov|NUM>: maximum coverage threshold. Will look for SAMPLE.Q${Q}.nodup.cov${maxcov}.coverage.txt [DEFAULT: 500]
4) <dir|CHR>:    folder with subdirectories for each sample [DEFAULT: current directory ]
```

**Depends on**
```
# R packages:
ggplot2
tidyr
```


**Example**
```
collect.coverage.stats.R samples.txt 20
```


## 4) filter.visual.coverages.sh

**Usage**
```
filter.visual.coverages.sh -s <file> -t <file> -r <file> -a <numeric fraction> -b <numeric fraction> -c <positive integer> \
                           -d <positive integer> -e <positive integer> -f <numeric fraction> -p <numeric fraction>
```

**Arguments**
```
# Required
-s         Path to samples file. Header and tab-separation expected.
          
           Sample IDs must be in the FIRST column. These must match (a subset of) sample names in the mapping 
           stats passed via -t.
          
           Group IDs can be specified in the SECOND column (if not specified, all 
           samples are assumed to constitute one group).
          
           The group ID is used to apply region filtering criteria 4-9 within all considered groups, to determine regions
           passing the filtering criteria in all groups.
          
           Samples that do not belong to any specified group (second column empty or 'NA') will be displayed in summary
           plots but will not be considerd during region filtering. 
          
           Additional columns are ignored.
          
          
-t         Path to coverage statistics. Header and tab-separation expected.

           Sample IDs must be in the FIRST column. Coverage statistics must be in the following columns as defined in
           filter.visual.coverages.R.
          
           Only coverage statistics of samples passed via -s will be used. A Warning or Stop is issued if there
           are mismatches.
          

-r         Path to region reference sequences. FASTA format expected. Used to correlate alignment stats with 
           reference sequence lengths and GC content.
          
           Only target regions passed via -t will be considered. A Warning or Stop is issued if there are mismatches.


# Optional [DEFAULT]
           # The first two filters take absolute thresholds and aim to remove poorly sequenced samples or regions:
- a  [0.3] minimum fraction of regions with at least one mapped read in a sample (filters samples)
- b  [0.3] minimum fraction of samples with at least one mapped read in a region (filters target regions)
 
           # The next four filters take thresholds...
- c  [500] minimum BWA-MEM alignment length
- d   [10] minimum average coverage in the aligned region
- e [1000] maximum average coverage in the aligned region
- f  [0.5] minimum alignment fraction (BWA-MEM alignment length divided by target region length)

           # ...that need to be met in a specified *fraction* of samples in each considered taxon group:
- p  [0.9]  minimum fraction of samples in each taxon group that need to pass filters c-f in order to keep a region

```

**Depends on**
```
# R packages:
ggplot2
ape
grid
VennDiagramm
tidyr
```


**Example**
```
filter.visual.coverages.sh -s samples.mapfile.txt -t coverage_stats.txt -r reference.fasta \
                           -a 0.2 -b 0.4 -c 1 -d 8 -e 1000 -f 0 -p 0.4
```

## Continue
[➜ Continue to Step 2](Step2_Sequence_Assembly.md)
