[← Back to Tutorial Home](../)

# Preprocessing Step 0.1: Sequence Quality Control

## 1) [run.fastqc.sh](https://github.com/scrameri/CaptureAl/wiki/run.fastqc.sh)

Runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for 12 samples given in [samples.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/samples.txt) in parallel, using sequences in [adapters.fasta](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/adapters.fasta) to estimate adaptor content.

**Example**
```
run.fastqc.sh -s samples.txt -d NovaSeq-run1_raw -a adapters.fasta -t 12
```

## 2) [plot.fastqc.R](https://github.com/scrameri/CaptureAl/wiki/plot.fastqc.R)

Visualize FastQC results as a PDF.

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) already provides a visualization of the results via .html files. However, if you work with many samples, you'd need to inspect dozens or hundreds of such files manually. CaptureAl provides a concise visualization of up to hundred or more `*_fastqc.zip` result files, in a single PDF file of 18 pages. Some plots show the distribution of key statistics across all samples, while other plots show detailled statistics on a per-sample basis using `ggplot2` and facets for a direct comparison across samples.

**Example**
```bash
cd NovaSeq-run1_raw
ls -1 fastqc/*_fastqc.zip > samples.fastqc.txt
plot.fastqc.R samples.fastqc.txt fastqc_raw.pdf 18 18
```

## Continue
[➜ Continue to Read Trimming](Step0.2_Read_Trimming.md)

