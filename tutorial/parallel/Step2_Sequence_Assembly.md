[← Back to Step 1](Step1_Read_Mapping.md)


# STEP 2

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step2.png)


## 1) [extract.readpairs.sh](https://github.com/scrameri/CaptureAl/wiki/extract.readpairs.sh)

This creates an output directory with a subdirectory for each sample in [samples.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/tutorial/data/samples.txt), containing read pairs extracted from `*fastq.gz` files in `NovaSeq-run1_trimmed`. Only read pairs where one or both reads mapped to the reference sequence of a locus in [loci.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/tutorial/data/loci.txt) are extracted, as determined by mapping results in `NovaSeq-run1_mapped` and `-Q`.

It is highly recommended to run this on a local scratch, due to the large number of files written.

**Example**
```
extract.readpairs.sh -s samples.txt -l loci.txt -d NovaSeq-run1_trimmed -m NovaSeq-run1_mapped \
                     -b NovaSeq-run1_mapped/SAMPLE/SAMPLE.bwa-mem.sorted.Q10.nodup.bam -Q 10 -t 20
```

## 2) [run.dipspades.sh](https://github.com/scrameri/CaptureAl/wiki/run.dipspades.sh)

This assembles the extracted reads of each locus and each sample using dipSPAdes.

It is highly recommended to run this on a local scratch, due to the large number of files written.

**Example**
```
run.dipspades.sh -s samples.txt -r NovaSeq-run1_extracted -t 20
```

## Continue
[➜ Continue to Step 3](Step3_Orthology_Assessment.md)
