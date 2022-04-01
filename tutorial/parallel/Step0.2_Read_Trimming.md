[← Back to Sequence Quality Control](Step0.1_Sequence_Quality_Control.md)

# Preprocessing Step 0.2: Read Trimming

## 1) [trim.fastq.sh](https://github.com/scrameri/CaptureAl/wiki/trim.fastq.sh)

[samples.fabaceae.12.txt]()

```
trim.fastq.sh -s samples.fabaceae.12.txt -a illumina.truseq.indexing.adaptors -r NovaSeq-run1_raw -x '_R1.fastq.gz,_R2.fastq.gz' -t 20
```

## Continue
[➜ Continue to Step 1](Step1_Read_Mapping.md)
