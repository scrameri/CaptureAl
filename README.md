# CaptureAl
Pipeline streamlining analysis of **target enrichment** sequencing data


# Tutorial using toy dataset [in prep.]

### Set Parameters
```bash
# input directories
scratch="/cluster/home/crameris/scratch"
raw="~/home/Dalbergia/uce/seq-rawdata/novaseq-run9-10_Dec.2021"
batch=$(basename $raw)                          
mappingdir="mapping/${batch}"

# input files
adapter="/cluster/work/gdc/people/crameris/Dalbergia/uce/seq-qualfiltered/adapters.txt"
ref="/cluster/work/gdc/people/crameris/Dalbergia/uce/references/consDalbergia_4c_2396.fasta"

# output directories
fastqc_raw="fastqc_raw"
fastqc_trim="fastqc_trim"

threads=10
```

### PREPROCESSING READS
#### Define samples to process

You should tell the pipeline which files with raw sequencing reads (typically ending in `*fastqc.gz` or `*fq.gz`) should be processed. Some analysis steps, such as *read trimming* and *mapping* can be done separately for each file, while later steps such as *alignment* require that a collection of samples are analyzed up to a certain step.

Let's start by going to our working directory `$scratch`. This directory will be used to write a lot of output, and scratch directories are typically the best-suited to do so efficiently. However, bear in mind that data is stored only temporarily on scratch directories, with no backup.

Let's also peek into the directory containing files with raw reads (`$raw`) and extract the *basename* of every sample located there. If the data is paired-end, a sample will have forward and reverse reads, and downstream programs such as *Trimmomatic* need to know which files belong together. Files with the same basename are assumed to represent forward and reverse reads of the same DNA library.

```bash
cd $scratch
find ${raw} -maxdepth 1 -regextype egrep -regex '.*[_.]R1[_.].*' |sed 's!.*/!!' |sed 's/[.][/]//' |sed 's/[_.]R1[_.].*//' |sort |uniq > samples.txt
```

Have a look at `samples.txt` to see if the expected basenames are there. There should be half as many lines in `samples.txt` as `*fastq.gz` files in `$raw`.


#### Run FastQC on raw reads

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a program to summarize and visualize key statistics from read data, such as the number of sequencing reads, read quality, and GC content. Once the program has been installed on your system, you should be able to run it using the `fastqc` command.

The most important arguments of `fastqc` are `-t` (number of threads), `-o` (output directory, needs to exist) and `-a` (non-standard FASTA file with adapter sequences).

Depending on your computing environment, there is the **BSUB solution** and the **GNU parallel** solution to run programs for many samples in parallel. The first depends on `bsub` and is suitable for execution on large high-performance clusters, while the second uses GNU `parallel` and is suitable for execution on clusters or personal computers. I'll provide both solutions where possible. Scripts for both solutions are available in the CaptureAl repository.

This will create an output directory `$fastqc_raw` and write two output files (`*_fastqc.html` and `*_fastqc.zip`) for each input file (`*.gz`).


**GNU Parallel solution**
```bash
if [ ! -d $fastqc_raw ] ; then mkdir $fastqc_raw ; fi
parallel -a <(ls -1 ${raw}/*fastq.gz) -j $threads fastqc -t 1 -o $fastqc_raw
```

**BSUB solution**
```bash
bsub < bsub.fastqc.sh
```

#### Visualize FastQC results as a PDF

*FastQC* already provides a visualization of the results via .html files. However, if you work with many samples, you'd need to inspect dozens or hundreds of such files manually. CaptureAl provides a concise visualization of up to hundred or more `*_fastqc.zip` result files, in a single PDF file of 18 pages. Some plots show the distribution of key statistics across all samples, while other plots show detailled statistics on a per-sample basis using `ggplot2` and facets for a direct comparison across samples.

The FastQC plotting function is `plot.fastqc.R`. It's arguments can be seen by executing the function without arguments. They are:

1) file with paths to `*_fastqc.zip` files to process
2) PDF name (default: "fastqc.pdf")
3) PDF height (default: 12 inches)
4) PDF width (default: 12 inches)

```bash
cd $fastqc_raw
ls -1 *_fastqc.zip > samples.fastqc.txt
plot.fastqc.R samples.fastqc.txt fastqc_raw.pdf 18 18 > /dev/null
cd ../
```

#### Trim raw reads using Trimmomatic

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a program to trim raw reads. Once the program has been installed on your system, you should be able to run it using the fastqc command.

Read trimming is the removal of (parts of) reads with low sequencing quality, small read length, or adapter contaminations. The program takes many arguments, these are the most important ones:

1) ILLUMINACLIP:${adapter}:2:30:10 
2) LEADING:3 
3) TRAILING:3 
4) SLIDINGWINDOW:4:15 
5) MINLEN:50

This will create an output directory `$batch` and write a `*.trim1.*.fastq.gz` (forward trimmed reads) and a `*.trim2.*.fastq.gz` (reverse trimmed reads) file for every sample in `samples.txt`. Unpaired reads (`*.U1.*fastq.gz`, `*.U2.*fastq.gz`) and log files (`*.log`, `*.err`) with process prints and errors will be located in the `logs` subdirectory.

***GNU Parallel solution***
```bash
if [ ! -d $batch ] ; then mkdir $batch ; fi
cd $batch
cp $scratch/samples.txt .
trim.fastq.sh -s samples.txt -a $adapter -r $raw -x '_R1.fastq.gz,_R2.fastq.gz' -t $threads
```

***BSUB solution***
```bash
bsub < bsub.trimmomatic.sh
```

#### Run FastQC on trimmed reads
Let's run *FastQC* on the trimmed reads now.

***GNU Parallel solution***
```bash
if [ ! -d $fastqc_trim ] ; then mkdir $fastqc_trim ; fi
parallel -a <(ls -1 *fastq.gz) -j $threads fastqc -t 1 -o $fastqc_trim
```

***BSUB solution***
```bash
bsub < bsub.fastqc.sh
```

Let's now visualize the *FastQC* results of the trimmed reads, and compare them to the results of the raw reads.

```bash
cd $fastqc_trim
ls -1 *_fastqc.zip > samples.fastqc.txt
plot.fastqc.R samples.fastqc.txt fastqc_trimmed.pdf 18 18 > /dev/null
```





### STEP 1: READ MAPPING
#### mapping
This runs ```bwa mem``` for samples in ```samples.txt``` using ```reference.fasta``` as reference sequences and 4 \* 5 = 20 threads in parallel. The input are trimmed reads located in the ```path/to/reads``` directory. Input is specified as *paired-end* reads, with two files per sample, ending in ```.trim1.fastq.gz``` and ```.trim2.fastq.gz```, respectively.
```
fix.fasta.headers.sh $ref
run.bwamem.sh -s samples.txt -r $ref -e .trim1.fastq.gz,.trim2.fastq.gz -T 10 -Q 20 -d path/to/reads/ -t 4
```

#### coverage analysis
This computes coverage statistics for samples in ```samples.txt``` (for reads mapped with quality 20), using 4 threads in parallel.
```
get.coverage.stats.sh -s samples.txt -Q 20 -t 4
```

This collects coverage statistics for samples in ```samples.txt``` (for reads mapped with quality 20), and writes them all to one file.

```
collect.coverage.stats.R samples.txt 20
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


