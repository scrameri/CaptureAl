[← Back to Step 4](Step4_Sample_and_Locus_Filtering.md) ................................................................................................................................................................. [➜ Continue to Step 6](Step6_Merge_Overlapping_Alignments.md)


# STEP 5

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step5.png)


## 1) create.multifastas.parallel.sh

**Usage**
```
create.multifastas.parallel.sh -s <file> -l <file> -d <directory> -t <positive integer>
```

**Arguments**
```
# Required
-s        Path to file with sample basenames (without header or '>').
-l        Path to file with locus names (without header or '>').
-d        Path to directory with exonerate results. This directory should contain a subdirectory for each sample,
          with FASTA files with the best-matching contigs.

# Optional [DEFAULT]
-t  [15]  Number of parallel threads.
```

**Details**
```
This script creates an output directory of the form `multifasta.${nind}.${nloc}`, where ${nind} is the number of 
samples and ${nloc} is the number of loci.
```

**Depends on**
```
--
```


**Example**
```
create.multifastas.parallel.sh -s samlples.txt -l loci.txt -d NovaSeq-run1_exonerate -t 20
```


## 2) align.multifastas.parallel.sh


**Usage**
```
align.multifastas.parallel.sh -d <directory> -m <string> -t <integer>
```

**Arguments**
```
# Required
-d                Path to directory with multifasta files. This directory should contain a FASTA file for each target region,
                  each with unaligned contigs of multiple samples.

# Optional [DEFAULT]
-m [`localpair`]  MAFFT alignment model. Choose between 'localpair', 'globalpair' or 'genafpair'. See MAFFT documentation for
                  further details.
-t         [4]    Number of parallel threads.
```

**Details**
```
This script creates an output directory of the form `mafft.${nind}.${nloc}`, where ${nind} is the number of 
samples and ${nloc} is the number of loci.
```


**Depends on**
```
mafft
```


**Example**
```
align.multifastas.parallel.sh -d multifasta.63.2396 -m 'localpair' -t 20
```


## 3) trim.alignment.ends.parallel.sh


**Usage**
```
trim.alignment.ends.parallel.sh -s <file> -d <directory> -c <numeric fraction> -n <numeric fraction> -t <positive integer> -v
```

**Arguments**
```
# Required
-s            File with sample names in FIRST column. Header and tab-separation expected. Any additional columns are ignored.
-d            Path to directory with raw alignments. This directory should contain a FASTA file for each target region,
              each with aligned contigs of multiple samples.

# Optional [DEFAULT]
-c    [0.5]   Completeness parameter. Alignments are trimmed at both ends until an alignment site has nucleotides in at least
              the specified fraction of aligned sequences. Both thresholds (-c and -n) need to be reached for trimming to stop.
-n   [0.25]   Maximum nucleotide diversity parameter (i.e., the sum of the number of base differences between sequence pairs,
              divided by the number of comparisons). Alignments are trimmed at both ends until an alignment site shows
              a nucleotide diversity of 0.25 or lower. Both thresholds (-c and -n) need to be reached for trimming to stop.
-m    ['-']   Gap character. This character is interpreted as missing data or a gap when using the -c and -n filters above.
-v  [false]   FLAG, if turned on, the alignment trimming will be visualized as a PDF (recommended for few alignments only).
-w     [15]   Width of output PDF file.
-h      [7]   Height of output PDF file.
-t      [4]   Number of parallel threads.

```

**Details**
```
This script creates an output directory of the form `<inputdirectory>.c${c}.d{$n}`, where ${c} is the completeness parameter
and ${n} is the nucleotide diversity parameter.
```

**Depends on**
```
# R package:
ape
```


**Example**
```
# no visualization
trim.alignment.ends.parallel.sh -s mapfile.dalbergia.txt -d mafft.63.2396 -c 0.5 -n 0.25 -t 20

# with visualization
trim.alignment.ends.parallel.sh -s mapfile.dalbergia.txt -d mafft.63.2396 -c 0.5 -n 0.25 -t 20 -v
```

## 4) trim.alignments.parallel.sh


**Usage**
```
trim.alignments.parallel.sh -s <file> -d <directory> -c <numeric fraction> -z <positive integer> \
                            -n <numeric fraction> -S <positive integer> -m <string> -t <positive integer> \
                            -w <positive integer> -h <positive integer> -i -v
```

**Arguments**
```
# Required
-s            File with sample names in FIRST column. Header and tab-separation expected. Any additional columns are ignored.
-d            Path to directory with raw or end-trimmed alignments. This directory should contain a FASTA file for each
              target region, each with aligned contigs of multiple samples.

# Optional [DEFAULT]
-c    [0.3]   Completeness parameter. Any alignment site with nucleotides in less than this specified fraction of aligned
              sequences is removed. Aligned sequences are defined as any row in the alignment with nucleotide data (samples
              where a locus is entirely missing (see -m option) are ignored).
-z     [20]   Window size parameter. Potential mis-assemblies or mis-alignments in each sequence are identified using a
              sliding window approach with this specified window size (in number of bases). 
-S      [1]   Step size parameter. Potential mis-assemblies or mis-alignments in each sequence are identified using a
              sliding window approach with this specified step size (in number of bases). 
-n    [0.5]   Conservation parameter. Entire windows are successively trimmed at contig ends if more than this fraction of
              the nucleotides in the conserved part of the window deviate from the alignment consensus. 
              A conserved part of each window is defined as the alignment sites with nucleotides in at least 20% of samples,
              and where the frequencies of minor alleles are all below 30% without considering gaps.
              By default, the sliding window approach stops if a successive window survives the trimming (see -i option).
-m    ['-']   Gap character. This character is interpreted as missing data or a gap when using the -c and -n filters above.
-i  [false]   FLAG, if turned on, the sliding window approach is not only applied to contig ends, but extended to internal
              regions of the alignment (no stopping criterion used).
-v  [false]   FLAG, if turned on, the alignment trimming will be visualized as a PDF (recommended for few alignments only).
-w     [15]   Width of output PDF file.
-h      [7]   Height of output PDF file.
-t      [4]   Number of parallel threads.

```

**Details**
```
This script creates an output directory of the form `<inputdirectory>.c${c}.n{$n}`, where ${c} is the completeness parameter
and ${n} is the conservation parameter.
```

**Depends on**
```
# R package:
ape
```


**Example**
```
# no visualization
trim.alignments.parallel.sh -s mapfile.dalbergia.txt -d mafft.63.2396.c0.5.n0.25 -c 0.4 -z 20 -S 1 -n 0.5 -t 20

# with visualization
trim.alignments.parallel.sh -s mapfile.dalbergia.txt -d mafft.63.2396.c0.5.n0.25 -c 0.4 -z 20 -S 1 -n 0.5 -t 20 -v

# with visualization and internal trimming
trim.alignments.parallel.sh -s mapfile.dalbergia.txt -d mafft.63.2396.c0.5.n0.25 -c 0.4 -z 20 -S 1 -n 0.5 -t 20 -iv
```

## Continue
[← Back to Step 4](Step4_Sample_and_Locus_Filtering.md) ................................................................................................................................................................. [➜ Continue to Step 6](Step6_Merge_Overlapping_Alignments.md)
