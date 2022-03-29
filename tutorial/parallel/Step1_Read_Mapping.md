Go back to [previous step]()

**Usage**
```
trim.fastq.sh -s <sample file> -a <adapter file> -r <directory> -x <string> -t <integer>
```

**Arguments**
```
-s    sample file containing one field (w/o header) with the sample base names
-a    adapter file in FASTA format containing adapter sequences
-r    path to raw reads: will look in this folder for <${sample}${suffix1}> raw reads [DEFAULT: current working directory ]
-x    suffixes for fwd and rev reads: comma-separated, e.g. '_R1.fastq.gz,_R2.fastq.gz' if the fwd raw read file is ${sample}_R1.fastq.gz and the rev raw reads file is in ${sample}_R2.fastq.gz
-o    path to output directory [DEFAULT: current working directory ]
-t    threads (should not exceed 30)

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
