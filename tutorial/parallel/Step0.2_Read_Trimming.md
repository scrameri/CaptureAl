[← Back to Sequence Quality Control](Step0.1_Sequence_Quality_Control.md)

# Preprocessing Step 0.2: Read Trimming

## 1) [trim.fastq.sh](https://github.com/scrameri/CaptureAl/wiki/trim.fastq.sh)

This performs read trimming for 12 samples [samples.fabaceae.12.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/samples.fabaceae.12.txt) with paired-end data (files ending in `_R1.fastq.gz` and `_R2.fastq.gz`) located ind the `NovaSeq-run1_raw` directory.

It uses adapter sequences from [adapters.fasta](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/adapters.fasta) and performs trimming in parallel using 12 threads.

```
trim.fastq.sh -s samples.fabaceae.12.txt -a adapters.fasta -r NovaSeq-run1_raw -x '_R1.fastq.gz,_R2.fastq.gz' -t 12
```

## Continue
[➜ Continue to Step 1](Step1_Read_Mapping.md)
