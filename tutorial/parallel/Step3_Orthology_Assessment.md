[← Back to Step 2](Step2_Sequence_Assembly.md)


# STEP 3

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step3.png)


## 1) [select.best.contigs.per.locus.sh](https://github.com/scrameri/CaptureAl/wiki/select.best.contigs.per.locus.sh)


**Example**
```
select.best.contigs.per.locus.sh -s samples.txt -l loci.txt -r reference.fasta -d NovaSeq-run1_assembly -t 20
```


## 2) [combine.contigs.parallel.sh](https://github.com/scrameri/CaptureAl/wiki/combine.contigs.parallel.sh)

This reads the assembled contigs per sample and individual (from the directory passed via `-d`) and the best-matching contig and exonerate alingment statistics (from the directory passed via `-e`), and combines non-overlapping contigs if they pass the requirements passed via `-a` (minimum target alingment length) and `-c` (minimum normalmized alingment score)

**Example**
```
combine.contigs.parallel.sh -s samples.txt -d NovaSeq-run1_assembly -e NovaSeq-run1_exonerate -a 80 -c 2 -t 20
```


## 3) [collect.exonerate.stats.R](https://github.com/scrameri/CaptureAl/wiki/collect.exonerate.stats.R)

This collects all relevant assembly statistics of all samples. These statistics are needed for the upcoming pipeline step 4 (sample and target region filtering). The file `loci_stats.txt` is the main requirement for step 4, the file `loci_contignumbers.txt` can be plotted (see 4) below).

**Example**
```
collect.exonerate.stats.R samples.txt NovaSeq-run1_exonerate
```


## 4) [plot.contig.numbers.R](https://github.com/scrameri/CaptureAl/wiki/plot.contig.numbers.R)

Visualize the file `loci_contignumber.txt` produced in the previous step. This should give an idea on how many target regions are represented by single vs. multiple contigs (capture specificity), and how many target regions have no contig (capture sensitivity).

If a taxon group map file such as [mapfile.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/mapfile.txt) is provided (second argument), the resulting box plot is plotted by taxon groups.

**Example**
```
plot.contig.numbers.R loci_contignumbers.txt mapfile.txt
```


## Continue
[➜ Continue to Step 4](Step4_Sample_and_Locus_Filtering.md)
