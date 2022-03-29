Go back to [previous step]()

## Overview
![Step1.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/Step1.png)


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

# Optional
-r    path to raw reads: will look in this folder for <${sample}${suffix1}> raw reads [DEFAULT: current working directory ]
-x    comma-separated suffixes for forward (fwd) and reverse (rev) reads, e.g. '_R1.fastq.gz,_R2.fastq.gz' if the fwd raw read file is ${sample}_R1.fastq.gz and the rev raw reads file is in ${sample}_R2.fastq.gz [DEFAULT:'_R1.fastq.gz,_R2.fastq.gz']
-o    path to output directory [DEFAULT: current working directory ]
-t    number of threads for parallel computations [DEFAULT: 16]
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

Go to [next step]()
