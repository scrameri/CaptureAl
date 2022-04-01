[← Back to Tutorial Home](../)

# Preprocessing Step 0.1: Sequence Quality Control

## 1) [1_fastqc_raw.sh](https://github.com/scrameri/CaptureAl/wiki/1_fastqc_raw.sh)

**Depends on**
```
adapters.txt
```


**Example**
```
cd NovaSeq-run1_raw
../1_fastqc_raw.sh
```

## 2) Visualize FastQC results as a PDF

**Usage**
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) already provides a visualization of the results via .html files. However, if you work with many samples, you'd need to inspect dozens or hundreds of such files manually. CaptureAl provides a concise visualization of up to hundred or more `*_fastqc.zip` result files, in a single PDF file of 18 pages. Some plots show the distribution of key statistics across all samples, while other plots show detailled statistics on a per-sample basis using `ggplot2` and facets for a direct comparison across samples.

The FastQC plotting function is [plot.fastqc.R](https://github.com/scrameri/CaptureAl/wiki/plot.fastqc.R). It's arguments can be seen by clicking on the link (to the wiki page) or by executing the function without arguments.


**Example**
```bash
ls -1 fastqc/*_fastqc.zip > samples.fastqc.txt
plot.fastqc.R samples.fastqc.txt fastqc_raw.pdf 18 18
```

## Continue
[➜ Continue to Read Trimming](Step0.2_Read_Trimming.md)

