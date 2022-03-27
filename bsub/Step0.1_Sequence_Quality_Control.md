# Step 0.1: Sequence Quality Control

### Define working directories

It is highly recommended to work on a scratch disk if you have access to one, as this will speed up analyses and prevent automatic backup of heavy unzipped or unimportant temporary files generated during the analysis. However, bear in mind that data is stored only temporarily on scratch, and should be moved to a directory with backup once the analysis is finished, e.g. using the `mv` or `rsync` commands.

This is how I set my working directories

```bash
# my scratch directory
scratch="/cluster/home/crameris/scratch"

# path to folder with raw reads (where the raw data is located)
raw="${scratch}/CaptureAl/seq-rawdata/tutorial"

# path to folder with quality-filtered reads (does not exist yet)
trimmed="${scratch}/CaptureAl/seq-qualfiltered/tutorial"

# path to reference .fasta file
ref="${scratch}/reference.fasta"

```

### Define samples to process

Steps 0-3 can be executed for any collection of samples separately (in parallel), while later steps 4-7 require that a collection of samples is analyzed up to a certain step, and parallelization is over genetic regions (loci). 

In most steps, you can tell the pipeline which samples or loci should be processed. By default, all samples in ${raw} and all loci in ${ref} are being processed.

Let's start by changing to our working directory `${raw}`. 

```bash
cd ${raw}
```

The *CaptureAl* pipeline needs to know which samples should be processed. If the sequence data is paired-end (as in most target capture experiments), a sample will have forward reads (file name ending in `_R1.fastq.gz`) and reverse reads (file name ending in `_R2.fastq.gz). Downstream analysis steps need to know which files belong to the same individual DNA library. These file pairs should have a common *base name* coding for the individual DNA library, such as the base name `SH0356` for files `SH0356_R1.fastq.gz` and `SH0356_R2.fastq.gz`. 

We can extract the *base name* of every individual DNA library located in our directory as follows:

```bash
find ${raw} -maxdepth 1 -regextype egrep -regex '.*[_.]R1[_.].*' |sed 's!.*/!!' |sed 's/[.][/]//' |sed 's/[_.]R1[_.].*//' |sort |uniq > samples.txt
```

Have a look at `samples.txt` and see if the expected sample base names are there. There should be half as many lines in `samples.txt` as `*fastq.gz` files in `${raw}`. You can check this by typing

```bash
wc -l samples.txt # gives number of lines (sample base names) in file <samples.txt>
ls -1 *.fastq.gz | wc -l # gives number of *.fastq.gz files in your working directory
```


### Run FastQC on raw reads

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a program to summarize and visualize key statistics from read data, such as the number of sequencing reads, read quality, and GC content. Once the program has been installed on your system, you should be able to run it using the `fastqc` command.

The most important arguments of `fastqc` are `-t` (number of threads), `-o` (output directory) and `-a` (optional FASTA file with non-standard adapter sequences).

This will create subdirectory `fastqc` and write two output files (`*_fastqc.html` and `*_fastqc.zip`) for each input file (`*.fastq.gz`).

```bash
bsub < bsub.fastqc.sh
```

### Visualize FastQC results as a PDF

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) already provides a visualization of the results via .html files. However, if you work with many samples, you'd need to inspect dozens or hundreds of such files manually. CaptureAl provides a concise visualization of up to hundred or more `*_fastqc.zip` result files, in a single PDF file of 18 pages. Some plots show the distribution of key statistics across all samples, while other plots show detailled statistics on a per-sample basis using `ggplot2` and facets for a direct comparison across samples.

The FastQC plotting function is `plot.fastqc.R`. It's arguments can be seen by executing the function without arguments. They are:

1) file with paths to `*_fastqc.zip` files to process
2) PDF name (default: "fastqc.pdf")
3) PDF height (default: 12 inches)
4) PDF width (default: 12 inches)

```bash
ls -1 fastqc/*_fastqc.zip > samples.fastqc.txt
plot.fastqc.R samples.fastqc.txt fastqc_raw.pdf 18 18
```

### Next steps
To get to the next steps, follow the [Read Trimming](https://github.com/scrameri/CaptureAl/blob/master/bsub/Step0.2_Read_Trimming.md) tutorial.
