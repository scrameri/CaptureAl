[← Back to Step 6](Step6_Merge_Overlapping_Alignments.md)


# STEP 7

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step7.png)


## 1) [get.group.consensus.sh](https://github.com/scrameri/CaptureAl/wiki/get.group.consensus.sh)

This creates a consensus sequence for every alignment in the directory passed via `-d`. Taxon groups specified in [mapfile.txt](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/mapfile.txt) are used to make the consensus sequence representative of all taxon groups, irrespective of their rarity.

The resulting consensus sequences can be used as target region reference sequences for another iteration of mapping, assembly, and alignment if needed.

**Example**
```
get.group.consensus.sh -s mapfile.txt -d mafft.63.2396.c0.5.d0.25.c0.4.n0.5 -m 1 -b 0.01 \
                       -z '.all.aln.etr.itr.cons' -t 20 -gnv
```

## Continue
[➜ Reiterate beginning with Step 1](Step1_Read_Mapping.md)
