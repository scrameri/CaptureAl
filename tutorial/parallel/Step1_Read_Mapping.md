[← Back to Read Trimming](Step0.2_Read_Trimming.md)


# STEP 1

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step1.png)


## 1) run.bwamem.sh

**Usage**
```
run.bwamem.sh -s <file> -r <file> -e <string> -d <directory> -T <positive integer> -Q <positive integer> -o <directory> -t <positive integer>
```

**Arguments**
```
# Required
-s        File with sample names (without header or '>')
-r        File with reference sequences (FASTA format)
-e        Sample file extension(s). E.g. '.fasta' or '.trim.fq.gz' [unpaired] or '.trim1.fastq,.trim2.fastq' [paired].
          The program then interprets if data is unpaired or paired. Separate file extensions of file pairs with a ','.

# Optional [DEFAULT]
-d  [pwd] Path to input reads.
-T  [10]  Minimum bwa mem alignment score, passed to -T parameter of bwa mem.
-Q  [20]  Minimum mapping quality, used to filter by the fifth field / MAPQ column in BAM files. Must be Q >= T.
-o  [pwd] Path to output directory. A folder will be created if it does not exist.
-t  [3]   Number of samples processed in parallel. Can be between 1 (uses ${cpu} CPU cores in total) and 6 (uses 6*${cpu} CPU cores in total, cpu=4).
```

**Details**
```
-s  The sample file must contain the sample basenames (i.e., sample name without file extensions specified in -e).
    Use a single line for paired reads in two files, e.g. the sample basename for files 'SH598_S16.trim1.fastq.gz'
    and 'SH598_S16.trim2.fastq.gz' would be 'SH598_S16' if -e is set to '.trim1.fastq.gz,.trim2.fastq.gz'.

-Q  A value above 0 will filter reads with multiple mappings. Only reads passing the -Q filter will be written to
    the output directory, after removing PCR duplicates using the picard MarkDuplicates tool.
```

**Depends on**
```
bwa
sambamba
java
python
picard-tools
r
```


**Example**
```
run.bwamem.sh -s samles.txt -r reference.fasta -e '.trim1.fastq.gz,.trim2.fastq.gz' -T 10 -Q 20 -d NovaSeq-run1_trimmed -t 4
```

## 2) get.coverage.stats.sh

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

## 3) collect.coverage.stats.R $s $Q

**Usage**
```
collect.coverage.stats.R <file> <positive integer>
```

**Arguments**
```
# Required
1) <sfile|CHR>:  path to sample file containing paths to sample directories with mapping results. No header expected. ;

# Optional
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
filter.visual.coverages.sh -s <file> -t <file> -r <file> -a <numeric fraction> -b <numeric fraction> -c <positive integer> -d <positive integer> -e <positive integer> -f <numeric fraction> -p <numeric fraction>
```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
filter.visual.coverages.sh -s samples.mapfile.txt -t coverage_stats.txt -r reference.fasta -a 0.2 -b 0.4 -c 1 -d 8 -e 1000 -f 0 -p 0.4
```

## Continue
[➜ Continue to Step 2](Step2_Sequence_Assembly.md)
