# CaptureAl
Pipeline streamlining analysis of **target enrichment sequencing** data


# Tutorial using toy dataset [in prep.]

### Check Installation
This repository provides scripts that use third-party software. These sofware tools must be executable from your computing environment:
- R and corresponding Rscript (tested on version 3.1.2, 3.6 and 4.0)
- python (tested on version 2.7.11)
- java (tested on version 1.8.0_73)
- perl (tested on version 5.16.3)
- fastqc (tested on version 0.11.4)
- trimmomatic (tested on version 0.35)
- bwa (tested on version 0.7.17)
- sambamba (tested on version 0.8.0)
- samtools (tested on version 1.2)
- bedtools (tested on version 2.28.0)
- picard-tools (tested on version 2.23.8)
- spades or dipspades (tested on version 3.14.1)
- exonerate (tested on version 2.4.0)
- freebayes (tested on version 1.3.4)

### Set Parameters
```bash
# input directories
scratch="/cluster/home/crameris/scratch"
raw="${scratch}/CaptureAl/seq-rawdata/tutorial"
trimmed="${scratch}/CaptureAl/seq-qualfiltered/tutorial"
ref="/cluster/home/crameris/reference.fasta"
mappingdir="${scratch}/mapping-reads-to-2396/${run}"
```

Depending on your computing environment, there is the **GNU parallel** solution and the **BSUB solution** to run programs for many samples or loci in parallel. The first uses GNU `parallel` and is suitable for execution on clusters or personal computers, and the second depends on `bsub` and is suitable for execution on large high-performance clusters. In each case, multiple so-called jobs are distributed and executed on different computing nodes.

The GNU parallel scripts take arguments via the available options, which are explained in the first section of each script (e.g. execute `somescript.sh -s samples.txt -t 5` to apply an analysis step to all samples specified in `samples.txt` using 5 threads). For reproducibility, automatically created log files will document which arguments were passed to each script.

The BSUB scripts have all their arguments set in the script (section `# arguments`), and they are therefore scripts and log files at the same time. Where appropriate, log files are still automatically created, which ensures reproducibility. The section `## Resource usage` is used to set the amount of computing nodes (threads), memory (in MB) and time (in hours) needed for each submitted job. These need to be set according to the amount of data analyzed. They should be as close as possible to the actual resources being used. It's always good practice to test the resources needed per job on a subset of e.g. 5 jobs before starting a bing sequence of jobs. This can be achieved by inspecting the lsf.o${job_ID} file, e.g. using the program `get.lsf.summary.sh lsf.o${job_ID}`, which displays a summary of average and maximum memory usage and computation time. This helps to set memory and time requirements for a big sequence of jobs, and prevents that submitted jobs are allocated too much resources, leading to inefficient use of shared computing power at the expense of other users, or too little resoruces, leading to premature job termination (automatic killing) with no results.

This tutorial uses BSUB scripts, but analogous scripts for both solutions are available in the CaptureAl repository.


### PREPROCESSING READS
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

#### Trim raw reads using Trimmomatic

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a program to trim raw reads. Once the program has been installed on your system, you should be able to run it using the fastqc command.

Read trimming is the removal of (parts of) reads with low sequencing quality, small read length, or adapter contaminations. The program takes many arguments, these are the most important ones:

1) ILLUMINACLIP:${adapter}:2:30:10 
2) LEADING:3 
3) TRAILING:3 
4) SLIDINGWINDOW:4:15 
5) MINLEN:50

This will create an output directory `${trimmed}` and write a `*.trim1.*.fastq.gz` (forward trimmed reads) and a `*.trim2.*.fastq.gz` (reverse trimmed reads) file for every sample in `samples.txt`. Unpaired reads (`*.U1.*fastq.gz`, `*.U2.*fastq.gz`) and log files (`*.log`, `*.err`) with process prints and errors will be located in the `logs` subdirectory.

```bash
bsub < bsub.trimmomatic.sh
```


#### Run FastQC on trimmed reads
Let's run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on the trimmed reads now. The `bsub.fastqc.sh` script is aware of your current directory and will produce a `fastqc` subfolder with results for the trimmed reads.

```bash
cd ${trimmed}
bsub < bsub.fastqc.sh
```

Let's now visualize the *FastQC* results of the trimmed reads, and compare them to the results of the raw reads.

```bash
ls -1 fastqc/*_fastqc.zip > samples.fastqc.txt
plot.fastqc.R samples.fastqc.txt fastqc_trimmed.pdf 18 18
```
You should see on page 10 that the per sequence quality scores have increased slightly, and on page 17 that the adapter content has decreased to very low levels now.


### STEP 1: READ MAPPING
#### mapping
This runs `bwa mem` for samples in `samples.txt` using `${ref}` as reference sequences. The input are trimmed reads located in the ${in} input directory. Input is specified as *paired-end* reads, with two files per sample, ending in `.trim1.fastq.gz` and `.trim2.fastq.gz`, respectively.

```
bsub < bsub.bwamem.sh
```
This creates an output directory ${mapping} with a subfolder for each sample, as well as a copy of the sample file `samples.txt`. In each sample subdirectory, you'll find the reference `*.fasta` file used to map against, the mapped reads in the `*.bam` files and corresponding `*.bam.bai` index files, as well as some basic mapping statistics ending in `*.flagstats.txt` and a `bwa.Q${Q}.log` log file.

#### coverage analysis
This computes coverage statistics for samples in `samples.txt` for reads mapped with quality ${Q} against all regions in `${ref}`.
```
cd ${mapping}
bsub < bsub.get.coverage.stats.sh
```

This collects mapping and coverage statistics for all samples in `samples.txt` using the mapping quality and coverage thresholds ${Q} and ${maxcov}, respectively.
```
collect.coverage.stats.sh -s samples.txt -Q ${Q} -m ${maxcov}
```

#### visualize and filter regions based on coverage statistics
This generates a heatmap and violin plots of the coverage analysis results. The visual output can be used to update the selection of adequate filtering thresholds. The ```mapfile.txt``` can be a file with a header and single column with sample base names, or it can contain a second column with group memberships. Group memberships are used to filter loci in each group separately, and to determine the overlap of remaining loci ```loci.txt```.
```
# Filtering thresholds
minploci=0.2 # min. proportion of regions recovered per sample (filters taxa)
minptaxa=0.4 # min. proportion of samples recovered per region (filters loci)
minlen=1     # minimum mapped length in .bam (filters loci)
mincov=8     # minimum average coverage in .bam (filters loci)
maxcov=1000  # maximum average coverage in .bam (filters loci)
minratio=0   # minimum target alignment fraction (alignment length / target length) (filters loci)
minfrac=0.4  # minimum fraction of samples conforming to the absolute locus filters (minlen, mincov, maxcov, minratio)

# Visualize coverage analysis and filter loci
filter.visual.coverages.R mapfile.txt coverage_stats.txt reference.fasta ${minploci} ${minptaxa} ${minlen} ${mincov} ${maxcov} ${minratio} ${minfrac}
````

### STEP 2: SEQUENCE ASSEMBLY
#### extract read pairs 
This extracts .fastq read pairs located in ```/path/to/reads``` for samples in ```samples.txt``` and loci in ```loci.txt```, and writes two separate .fastq files (one with forward reads, one with reverse reads) for each locus. One or both of the extracted reads per read pair mapped with mapping quality 20.

NOTE: This script should be executed on a local scratch, since many files are read and written, which is why the ```-m``` option can be used to redirect back to the working directory used in step 1.
```
extract.readpairs.sh -s samples.txt -r loci.txt -d path/to/reads/ -m path/to/mapping/ -Q 20 -t 4
```

This assembles reads in ```/path/to/extracted-reads/``` separately for each locus in ```loci.txt``` and each sample in ```samples.txt```, using 4 threads in parallel.

NOTE: This script should be executed on a local scratch, since many files are read and written, which is why the ```-r``` option can be used to redirect back to the working directory used in step 2.
```
run.dipspades.sh -s samples.txt -r /path/to/extracted-reads/ -t 4
```

### STEP 3: ORTHOLOGY ASSESSMENT

### STEP 4: SAMPLE AND LOCUS FILTERING

### STEP 5: TARGET LOCUS ALIGNMENT AND ALIGNMENT TRIMMING

### STEP 6: MERGE OVERLAPPING ALIGNMENTS

### STEP 7: CREATE REPRESENTATIVE REFERENCE SEQUENCES


