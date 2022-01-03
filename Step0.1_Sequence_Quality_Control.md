# Sequence Quality Control

#### Define samples to process

Some analysis steps, such as *read trimming* and *mapping* can be done separately for each sample with sequence data, while later steps such as *alignment* require that a collection of samples is analyzed up to a certain step. 

In most steps, you can tell the pipeline which samples or loci should be processed.

Let's start by changing to our working directory `${raw}`. It is highly recommended to work on a scratch disk, as this will speed up analyses and prevent automatic backup of heavy unzipped or unimportant temporary files generated during the analysis. However, bear in mind that data is stored only temporarily on scratch, and should be moved to a directory with backup once the analysis is finished, e.g. using the `mv` or `rsync` commands.

Since the pipeline needs to know which samples should be processed, we start by extracting the *basename* of every sample located in our working directory. If the data is paired-end, a sample will have forward (file name including `_R1_` or `_R1.`) and reverse (file name including `_R2_` or `_R2.`) reads, and these file pairs should have a common *basename*. Downstream programs such as *Trimmomatic* need to know which files belong to the same DNA library.

```bash
cd ${raw}
find ${raw} -maxdepth 1 -regextype egrep -regex '.*[_.]R1[_.].*' |sed 's!.*/!!' |sed 's/[.][/]//' |sed 's/[_.]R1[_.].*//' |sort |uniq > samples.txt
```

Have a look at `samples.txt` and see if the expected sample basenames are there. There should be half as many lines in `samples.txt` as `*fastq.gz` files in `${raw}`.


#### Run FastQC on raw reads

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a program to summarize and visualize key statistics from read data, such as the number of sequencing reads, read quality, and GC content. Once the program has been installed on your system, you should be able to run it using the `fastqc` command.

The most important arguments of `fastqc` are `-t` (number of threads), `-o` (output directory, needs to exist) and `-a` (non-standard FASTA file with adapter sequences).

This will create an output directory `${fastqc_raw}` and write two output files (`*_fastqc.html` and `*_fastqc.zip`) for each input file (`*.gz`).

```bash
cd ${raw}
bsub < bsub.fastqc.sh
```

#### Visualize FastQC results as a PDF

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
