# CaptureAl
Pipeline streamlining analysis of **target enrichment** sequencing data


# Tutorial using toy dataset

### STEP 1: READ MAPPING
#### mapping
This runs ```bwa mem``` for samples in ```samples.txt``` using ```reference.fasta``` as reference sequences and 4 \* 5 = 20 threads in parallel. The input are trimmed reads located in the ```path/to/reads``` directory. Input is specified as *paired-end* reads, with two files per sample, ending in ```.trim1.fastq.gz``` and ```.trim2.fastq.gz```, respectively.
```
run.bwamem.sh -s samples.txt -r reference.fasta -e .trim1.fastq.gz,.trim2.fastq.gz -T 10 -Q 20 -d path/to/reads/ -t 4
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

NOTE: This script should ideally be executed on a local scratch, since many files are read and written, which is why the ```-m``` option can be used to redirect back to the working directory used in step 1.
```
extract.readpairs.sh -s samples.txt -r loci.txt -d path/to/reads/ -m path/to/mapping/ -Q 20 -t 4
```

This assembles reads in ```/path/to/extracted-reads/``` separately for each locus in ```loci.txt``` and each sample in ```samples.txt```, using 4 threads in parallel.
```
run.dipspades.sh -s samples.txt -r /path/to/extracted-reads/ -t 4
```

### STEP 3: ORTHOLOGY ASSESSMENT

### STEP 4: SAMPLE AND LOCUS FILTERING

### STEP 5: TARGET LOCUS ALIGNMENT AND ALIGNMENT TRIMMING

### STEP 6: MERGE OVERLAPPING ALIGNMENTS

### STEP 7: CREATE REPRESENTATIVE REFERENCE SEQUENCES


