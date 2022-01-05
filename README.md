# CaptureAl Analysis Pipeline
Pipeline streamlining analysis of **target enrichment sequencing** data

### Download Capture Al scripts
To understand use this pipeline, clone the [CaptureAl repository](https://github.com/scrameri/CaptureAl) and follow the tutorial provided below. For each step, there is a limited number of scripts that have to be executed, and these scripts often execute other scripts located in the repository.

### Check Installation
This repository provides scripts (see [src](https://github.com/scrameri/CaptureAl/tree/master/src) directory) that use third-party software. These sofware tools must be executable from your computing environment:
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

### Parallel Computing
Depending on your computing environment, there is the **GNU parallel** solution and the **BSUB solution** to run the analysis steps for many samples or loci in parallel. The first uses the GNU command `parallel` and is suitable for execution on clusters or personal computers, and the second uses the `bsub` command and is suitable for execution on large high-performance clusters. In each case, multiple so-called jobs are distributed and executed on different computing nodes.

The GNU parallel scripts take arguments via the available options, which are explained in the first section of each script (e.g. execute `somescript.sh -s samples.txt -t 5` to apply an analysis step to all samples specified in `samples.txt` using 5 threads). For **reproducibility**, automatically created log files will document which arguments were passed to each script.

The BSUB scripts have all their arguments set in the script (section `# arguments`), and they are therefore scripts and log files at the same time. Where appropriate, log files are still automatically created, which ensures **reproducibility**. The section `## Resource usage` is used to set the amount of computing nodes (threads), memory (in MB) and time (in hours) needed for each submitted job. These need to be set according to the amount of data analyzed. They should be as close as possible to the actual resources being used. It's always good practice to test the resources needed per job on a subset of e.g. 5 jobs before starting a bing sequence of jobs. This can be achieved by inspecting the lsf.o${job_ID} file, e.g. using the program `get.lsf.summary.sh lsf.o${job_ID}`, which displays a summary of average and maximum memory usage and computation time. This helps to set memory and time requirements for a big sequence of jobs, and prevents that submitted jobs are allocated too much resources, leading to inefficient use of shared computing power at the expense of other users, or too little resoruces, leading to premature job termination (automatic killing) with no results.

This tutorial is based on **BSUB** submission scripts, but analogous scripts for both solutions are available in the CaptureAl repository, and extending this turorial to comprise a full example using both solutions is work in progress. To follow the tutorial step by step, view the numbered `.md` files, starting with [Sequence Quality Control](https://github.com/scrameri/CaptureAl/blob/master/Step0.1_Sequence_Quality_Control.md).
