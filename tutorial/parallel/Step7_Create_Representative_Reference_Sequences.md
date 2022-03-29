[← Back to Step 6](Step6_Merge_Overlapping_Alignments.md)


# STEP 7

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step7.png)


## 1) get.group.consensus.sh

**Usage**
```
get.group.consensus.sh -s <file> -d <directory> -m <numeric fraction> -b <numeric fraction> \
                       -a <string> -z <string> -t <positive integer> -gnv
```

**Arguments**
```
# Required
-s                  File with names of samples in FIRST column and taxon group names in SECOND column.
                    These samples / taxon groups are considered during consensus calculation.
                    Header and tab separation expected. Additional columns are ignored.
-d                  Path to directory with alignments.

# Optional
-m   [0.05]         Minimum allele frequency to call a IUPAC ambiguity. Interpreted as 'minimum allele count' if >1.
-b    [0.5]         Minimum base frequency to return a consensus (major allele or IUPAC ambiguity) instead of a gap.
                    Interpreted as 'minimum base count' if >1.
-g  [false]         FLAG, if turned on, then '-' characters will be ignored during consensus calculation.
-n  [false]         FLAG, if turned on, then 'N' characters will be ignored during consensus calculation.
-a  [false]         Prefix in consensus sequence name. If no prefix is desired, use 'FALSE' instead of ''.
-z  [false]         Suffix in consensus sequence name. If no suffix is desired, use 'FALSE' instead of ''.
-v  [false]         FLAG, if turned on, the alignment consensus will be visualized as a PDF
                    (recommended for few alignments only).
-w     [15]         Width of output PDF file.
-h      [7]         Height of output PDF file.
-t      [4]         Number of parallel threads.
```

**Depends on**
```
get.consensus.from.alignment.parallel.sh
create.multiconsensus.parallel.sh
align.multifastas.parallel.sh
```


**Example**
```
get.group.consensus.sh -s mapfile.dalbergia.txt -d mafft.63.2396.c0.5.d0.25.c0.4.n0.5 -m 1 -b 0.01 \
                       -z '.all.aln.etr.itr.cons' -t 20 -gnv
```

## Continue
[➜ Reiterate beginning with Step 1](Step1_Read_Mapping.md)
