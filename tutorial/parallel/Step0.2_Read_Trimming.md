[← Back to Sequence Quality Control](Step0.1_Sequence_Quality_Control.md)

# Preprocessing Step 0.2: Read Trimming

## 1) trim.fastq.sh

**Usage**
```
trim.fastq.sh -s <sample file> -a <adapter file> -r <directory> -x <string> -t <integer>
```

**Arguments**
```
# Required
-s    sample file containing one field (w/o header) with the sample base names
-a    adapter file in FASTA format containing adapter sequences

# Optional [DEFAULT]
-r  [pwd]     Path to raw reads folder. Will look in this folder for <${sample}${suffix1}> to find raw reads
-x  [details] Comma-separated suffixes for forward (fwd) and reverse (rev) reads.
              E.g. '_R1.fastq.gz,_R2.fastq.gz' if the fwd raw read file is ${sample}_R1.fastq.gz and the
              rev raw reads file is in ${sample}_R2.fastq.gz [DEFAULT:'_R1.fastq.gz,_R2.fastq.gz'].
-o  [pwd]     Path to output directory.
-t  [16]      Number of threads for parallel computations.
```

**Depends on**
```
java
trimmomatic
```


**Example**
```
trim.fastq.sh -s samples.fabaceae.12.txt -a illumina.truseq.indexing.adaptors -r NovaSeq-run1_raw -x '_R1.fastq.gz,_R2.fastq.gz' -t 20
```

## Continue
[➜ Continue to Step 1](Step1_Read_Mapping.md)
