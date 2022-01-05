# Step 2: Sequence Assembly

### Extract read pairs 
This extracts .fastq read pairs located in ```/path/to/reads``` for samples in ```samples.txt``` and loci in ```loci.txt```, and writes two separate .fastq files (one with forward reads, one with reverse reads) for each locus. One or both of the extracted reads per read pair mapped with mapping quality 20.

NOTE: This script should be executed on a local scratch, since many files are read and written, which is why the ```-m``` option can be used to redirect back to the working directory used in step 1.
```
extract.readpairs.sh -s samples.txt -r loci.txt -d path/to/reads/ -m path/to/mapping/ -Q 20 -t 4
```

### Assemble extracted reads
This assembles reads in ```/path/to/extracted-reads/``` separately for each locus in ```loci.txt``` and each sample in ```samples.txt```, using 4 threads in parallel.

NOTE: This script should be executed on a local scratch, since many files are read and written, which is why the ```-r``` option can be used to redirect back to the working directory used in step 2.
```
run.dipspades.sh -s samples.txt -r /path/to/extracted-reads/ -t 4
```

### Next steps
To get to the next steps, follow the [Orthology Assessment](https://github.com/scrameri/CaptureAl/blob/master/Step3_Orthology_Assessment) tutorial.
