[← Back to Read Trimming](Step0.2_Read_Trimming.md)


# STEP 1

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step1.png)


## 1) run.bwamem.sh

**Usage**
```
run.bwamem.sh -s <file> -r <file> -e <string> -T <positive integer> -Q <positive integer> -d <directory> -o <directory> -t <positive integer>
```

**Arguments**
```
# Required
-s  sample file
-r  reference fasta file
-e  sample file(s) extension(s). E.g. '.fasta' or '.trim.fq.gz' [unpaired] or '.trim1.fastq,.trim2.fastq' [paired]. The program then interprets if data is unpaired or paired. Separate file extensions of file pairs with a ','.

# Optional [DEFAULT]
-T  [10]  bwa mem -T alignment score [mapped reads will go to $TMPDIR]
-Q  [20]  minimum mapping quality (fifth field / MAPQ in .sorted.bam file) [high-quality reads will go to $TMPDIR]. Only deduplicated reads will go to ${out}/${name}/*bwa-mem.sorted.Q${Q}.nodup.bam. Must be Q >= T. If above 0, will filter reads with multiple mappings.
-d  [pwd] path to input reads
-o  [pwd] path to output directory (will be created if it does not exist)
-t  [3]   number of samples processed in parallel. Can be between 1 (uses ${cpu} CPU cores in total) and 6 (uses 6*${cpu} CPU cores in total, cpu=4 by default)
```

**Details**
```
The sample file (-s option) must contain the basenames of all samples to be mapped (basenames = filenames without file extensions, e.g. 'SH598_S16')

The sample file extensions (-e option) of files 'SH598_S16.trim1.fastq.gz' and 'SH598_S16.trim2.fastq.gz' might be '.trim1.fastq.gz,.trim2.fastq.gz'
-> only list one basename per sample for paired data (where the two files have the same basename, but different file extensions, see -e option)

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
run.bwamem.sh -s samles.txt -r reference.fasta -e '.trim1.fastq.gz,.trim2.fastq.gz' -T 10 -Q 10 -d NovaSeq-run1_trimmed -t 4
```

get.coverage.stats.sh -s $s -Q $Q -t 20
collect.coverage.stats.R $s $Q

samples.mapfile.txt

## Continue
[➜ Continue to Step 2](Step2_Sequence_Assembly.md)
