# CaptureAl
Pipeline streamlining analysis of **target enrichment** sequencing data


# Tutorial using toy dataset

### STEP 1: READ MAPPING
#### mapping
This runs ```bwa mem``` for samples in ```samples.txt``` using ```reference.fasta``` as reference sequences and 4 \* 5 = 20 threads in parallel. The input are trimmed reads located in the ```path/to/reads``` directory. It's specified as *paired-end* reads, with two files per sample, ending in ```.trim1.fastq.gz``` and ```.trim2.fastq.gz```, respectively.
```
run.bwamem.sh -s samples.txt -r reference.fasta -e .trim1.fastq.gz,.trim2.fastq.gz -T 10 -Q 20 -d path/to/reads/ -t 4
```
#### coverage analysis
This computes coverage statistics for samples in ```samples.txt``` mapped with quality 20, using 4 threads in parallel.
```
get.coverage.stats.sh -s samples.txt -Q 20 -t 4
```

### STEP 2: SEQUENCE ASSEMBLY

### STEP 3: ORTHOLOGY ASSESSMENT

### STEP 4: SAMPLE AND LOCUS FILTERING

### STEP 5: TARGET LOCUS ALIGNMENT AND ALIGNMENT TRIMMING

### STEP 6: MERGE OVERLAPPING ALIGNMENTS

### STEP 7: CREATE REPRESENTATIVE REFERENCE SEQUENCES


