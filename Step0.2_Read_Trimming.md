# Step 0.2: Read Trimming

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

### Next steps
To get to the next steps, follow the [Read Mapping](https://github.com/scrameri/CaptureAl/blob/master/Step1_Read_Mapping.md) tutorial.
