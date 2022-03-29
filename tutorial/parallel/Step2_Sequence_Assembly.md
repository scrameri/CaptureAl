[← Back to Step 1](Step1_Read_Mapping.md)


# STEP 2

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step2.png)


## 1) extract.readpairs.sh

**Usage**
```
extract.readpairs.sh -s <sample file> -l <locus file> -d <directory> -m <directory> -Q <integer> -t <integer>
```

**Arguments**
```
# Required
-s sample file
-l locus file
-d  []  absolute path to folder with quality-filtered reads
-m  []  absolute path to folder with mapping dirs

# Optional

-o  [seq-extracted] output directory (created if inexistent)
-Q  [10]            minimum mapping quality, as used for mapping using run.bwamem.sh
-b  [see details]   regex-path to BAM file. Use SAMPLE as wildcard.
-t  [2]             number of threads used
```

**Details**
```
-b    <SAMPLE> can be part of the string and will be replaced by the actual sample using regex,
      DEFAULT: <${mapdir}/SAMPLE/SAMPLE.bwa-mem.sorted.Q10.nodup.bam>."

```

**Depends on**
```
extract-reads-from-fastq.pl
```


**Example**
```
extract.readpairs.sh -s samples.txt -l loci.txt -d NovaSeq-run1_trimmed -m NovaSeq-run1_mapped \
                     -b NovaSeq-run1_mapped/SAMPLE/SAMPLE.bwa-mem.sorted.Q10.nodup.bam -Q 10 -t 20
```

## 2) run.dipspades.sh

**Usage**
```
run.dipspades.sh -s <sample file> -r <directory> -t <integer>
```

**Arguments**
```
# Required
-s sample file (will look for extracted reads inside ${name}.targets)
-r absolute path to extracted read pairs

# Optional
-t number of threads [DEFAULT: 15]: the parallelization is over reference sequences, not over individuals (these are assembled one after the other)
```

**Depends on**
```
submit-commands-var.pl
```


**Example**
```
run.dipspades.sh -s samples.txt -r NovaSeq-run1_extracted -t 20
```

## Continue
[➜ Continue to Step 3](Step3_Orthology_Assessment.md)
