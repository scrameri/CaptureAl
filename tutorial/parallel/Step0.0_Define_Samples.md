[← Back to Tutorial Home](../)

# Define samples to process
In multiple steps of the *CaptureAl* pipeline, it needs to know which samples should be processed. Steps 0-3 can be executed for any collection of samples separately (in parallel), while later steps 4-7 require that a collection of samples is analyzed up to [Step 3](Step3_Orthology_Assessment.md), and parallelization is over target regions (loci). 

If the sequence data is paired-end (as in most target capture experiments), a sample will have forward reads (file name ending in `_R1.fastq.gz`) and reverse reads (file name ending in `_R2.fastq.gz). Downstream analysis steps need to know which files belong to the same individual DNA library. These file pairs should have a common *base name* coding for the individual DNA library, such as the base name `SH0356` for files `SH0356_R1.fastq.gz` and `SH0356_R2.fastq.gz`. 

We can extract the *base name* of every individual DNA library located in our directory as follows:

```bash
raw=/path/to/raw-reads # use your path here

find ${raw} -maxdepth 1 -regextype egrep -regex '.*[_.]R1[_.].*' |sed 's!.*/!!' |sed 's/[.][/]//' |sed 's/[_.]R1[_.].*//' |sort |uniq > samples.txt
```
Depending on your file naming, you may need to change the regex above.

Have a look at `samples.txt` and see if the expected sample base names are there. There should be half as many lines in `samples.txt` as `*fastq.gz` files. You can check this by typing

```bash
wc -l samples.txt # gives number of lines (sample base names) in file <samples.txt>
ls -1 *.fastq.gz | wc -l # gives number of *.fastq.gz files in your working directory
```

## Continue
[➜ Continue to Read Trimming](Step0.2_Read_Trimming.md)
