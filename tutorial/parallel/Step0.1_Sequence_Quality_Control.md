[← Back to Tutorial Home](../)

# Preprocessing STEP 1: Sequence Quality Control

## 1) run.fastqc.sh

```
1_fastqc_raw.sh
```

**Arguments**
```
# Required
None. The script looks for all files ending in *fastq.gz and uses 10 parallel threads to run FASTQC using the adaptor file adapters.txt
```

**Depends on**
```
make.FASTQC.html.page.sh
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

The FastQC plotting function is `plot.fastqc.R`. It's arguments can be seen by executing the function without arguments. They are:

```
plot.fastqc.R <file> <outputfile> <integer> <integer>
```

**Arguments**
```
1) file with paths to `*_fastqc.zip` files to process
2) PDF name (default: "fastqc.pdf")
3) PDF height (default: 12 inches)
4) PDF width (default: 12 inches)
```

**Example**
```bash
ls -1 fastqc/*_fastqc.zip > samples.fastqc.txt
plot.fastqc.R samples.fastqc.txt fastqc_raw.pdf 18 18
```

## Continue
[➜ Continue to Read Trimming](Step0.2_Read_Trimming.md)

