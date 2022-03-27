# Step 1: Read Mapping
#### mapping

Set an environmental variable containing the directory name where raw reads are located. The same name can then be re-used to create matching output directories.

```bash
run=$(basename ${raw}) # name of directory with raw reads: tutorial in this example
ref=reference.fasta
Q=10
maxcov=100
```


This runs `bwa mem` for samples in `samples.txt` using `${ref}` as reference sequences. The input are trimmed reads located in the ${in} input directory. Input is specified as *paired-end* reads, with two files per sample, ending in `.trim1.fastq.gz` and `.trim2.fastq.gz`, respectively.


```
bsub < bsub.bwamem.sh
```
This creates an output directory ${mapping} with a subfolder for each sample, as well as a copy of the sample file `samples.txt`. In each sample subdirectory, you'll find the reference `*.fasta` file used to map against, the mapped reads in the `*.bam` files and corresponding `*.bam.bai` index files, as well as some basic mapping statistics ending in `*.flagstats.txt` and a `bwa.Q${Q}.log` log file.

#### coverage analysis
This computes coverage statistics for samples in `samples.txt` for reads mapped with quality ${Q} against all regions in `${ref}`.
```
mapping="${scratch}/mapping-reads-to-2396/${run}"
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
minptaxa=0.7 # min. proportion of samples recovered per region (filters loci)
minlen=1     # minimum mapped length in .bam (filters loci)
mincov=8     # minimum average coverage in .bam (filters loci)
maxcov=1000  # maximum average coverage in .bam (filters loci)
minratio=0   # minimum target alignment fraction (alignment length / target length) (filters loci)
minfrac=0.7  # minimum fraction of samples conforming to the absolute locus filters (minlen, mincov, maxcov, minratio)

# Visualize coverage analysis and filter loci
filter.visual.coverages.R mapfile.txt coverage_stats.Q${Q}.txt ${ref} ${minploci} ${minptaxa} ${minlen} ${mincov} ${maxcov} ${minratio} ${minfrac}
````

### Next steps
To get to the next steps, follow the [Read Assembly](https://github.com/scrameri/CaptureAl/blob/master/tutorial/bsub/Step2_Sequence_Assembly.md) tutorial.
