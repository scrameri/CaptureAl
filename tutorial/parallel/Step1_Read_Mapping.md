[← Back to Read Trimming](Step0.2_Read_Trimming.md)


# STEP 1

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step1.png)


## 1) [run.bwamem.sh](https://github.com/scrameri/CaptureAl/wiki/run.bwamem.sh)

This maps the reads of samples specified in [samples.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/samples.txt) located in the `NovaSeq_run1_trimmed` directory against a [reference.fasta](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/reference.fasta) FASTA file.

In this example, it assumes *paired-end* reads located in separate files, e.g. `BEN001.trim1.fastq.gz` and `BEN001.trim2.fastq.gz` for sample *BEN001*. It further processes 4 samples in parallel, but performs hyper-threading for each sample, using 5 times as many threads as specified in `-t`.

After mapping using a `-T` quality threshold, the BAM files are further filtered for high-quality mappings (`-Q` option), and PCR duplicates are removed using [Picard Tools](https://broadinstitute.github.io/picard/).

**Example**
```
run.bwamem.sh -s samples.txt -r reference.fasta -e '.trim1.fastq.gz,.trim2.fastq.gz' -T 10 -Q 20 \
              -d NovaSeq_run1_trimmed -t 4
```

## 2) [get.coverage.stats.sh](https://github.com/scrameri/CaptureAl/wiki/get.coverage.stats.sh)


This performs several consecutive steps, for multiple samples specified in [samples.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/samples.txt) in parallel. It 1) collects basic mapping statistics from `*flagstat` files, then 2) calculates covered lengths and average coverage for each reference sequence, then 3) creates a summary table with the number of regions exceeding 1, 10, 20, etc. fold average coverage in the aligned portion, and 4) writes all mapping and coverage statistics to two large tables..

**Example**
```
get.coverage.stats.sh -s samples.txt -Q 20 -t 20
```

## 3) [filter.visual.coverages.sh](https://github.com/scrameri/CaptureAl/wiki/filter.visual.coverages.sh)

This filters samples and loci (target regions) with poor sequence data, taking taxon groups into account and using filtering thresholds informed by comprehensive visualizations.

**Example**
```
filter.visual.coverages.sh -s samples.mapfile.txt -t coverage_stats.txt -r reference.fasta \
                           -a 0.2 -b 0.4 -c 1 -d 8 -e 1000 -f 0 -p 0.4
```

## Continue
[➜ Continue to Step 2](Step2_Sequence_Assembly.md)
