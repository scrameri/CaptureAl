[← Back to Step 5](Step5_Alignment_and_Alignment_Trimming.md)


# STEP 6

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step6.png)


## 1) [get.consensus.from.alignment.parallel.sh](https://github.com/scrameri/CaptureAl/wiki/get.consensus.from.alignment.parallel.sh)

This creates an output directory with a FASTA consensus sequence for each alignment in `mafft.63.2396.c0.5.d0.25`, considering samples in [samples.txt](raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/data/samples.txt).

**Example**
```
# Consensus with IUPAC ambiguity codes (m = 0.2), gaps, Ns, but no visualization
get.consensus.from.alignment.parallel.sh -s samples.txt -d mafft.63.2396.c0.5.d0.25 -m 0.2 -b 0.01 -t 20
```

If the `-g` flag is turned on, any site where the computed consensus is a gap (`-`) is removed from the consensus sequence (recommended to turn this on depending on the alignment software used to map against the consensus sequence). 

If the `-n` flag is turned on, any site where the computed consensus is `N` is removed from the consensus sequence (it's generally not recommended to turn this on).

If the `-v` flag is turned on, the alignment and consensus sequence are visualized as a PDF file (only recommended for a small alignment number).

```
# Consensus of the major alleles (m = 1), no gaps, no Ns, but with visualization
get.consensus.from.alignment.parallel.sh -s samples.txt -d mafft.63.2396.c0.5.d0.25 -m 1 -b 0.01 -t 20 -gnv
```


## 2) [rename.fasta.headers.R](https://github.com/scrameri/CaptureAl/wiki/rename.fasta.headers.R)

This fixes the FASTA sequence identifiers, which got altered with suffixes by mapping and trimming. 

Be careful, the renaming is potentially irreversible. Therefore, it is recommended to test it on a copy of your FASTA first.

**Example**
```
# Remove the suffix '.all.aln.etr.itr' from fasta sequence names
rename.fasta.headers.R consensus.fasta '.all.aln.etr.itr' FALSE FALSE

# Remove any text after '__' from fasta sequence names
rename.fasta.headers.R consensus.fasta FALSE FALSE '__'

# Remove any instance of '_LG_' or '_Scaffold_' and all text after '__' from fasta sequence names
rename.fasta.headers.R consensus.fasta '_LG_' '_Scaffold_' '__'

```


## 3) blast.vs.self.sh

**Usage**
```
blast.vs.self.sh <file>
```

**Arguments**
```
# Required
1)  Path to FASTA file.

# Optional
--
```

**Details**
```
A BLAST+ database is created from the FASTA file, and the same sequences are used as query and target in a
BLAST+ search with expect value (evalue) 1E-04.

The results are then filtered to retain only hits between pairs of DIFFERENT sequences,
and only ONE out of two reciprocal hits.
```

**Depends on**
```
 blast.fasta.seqs.sh
 filter.blast.vs.self.R
```


**Example**
```
blast.vs.self.sh consensus.fasta
```

## 4) find.overlapping.alignments.R

**Usage**
```
find.overlapping.alignments.R <file> <BOOLEAN> <string> <string>
```

**Arguments**
```
# Required
1) <bf|CHR>:                path to filtered blast vs. self results (*.blast.filtered) 
     
# Optional [DEFAULT] (if one or more are given, they need to be given in this order)
2) <check.lg|BOOLEAN>:      if TRUE, checks linkage group (LG) conformity based on query and subject identifiers, string.lg.before and string.lg.after [DEFAULT: TRUE]
3) <string.lg.before|CHR>:  unique string in locus name that comes just BEFORE the linkage group ID [DEFAULT: 'LG_']
4) <string.lg.after|CHR>:   string in locus name that comes just AFTER the linkage group ID [DEFAULT: '_']
5) <check.overlap|BOOLEAN>: if TRUE, checks whether hits are at alignment ends based on q/s.start/end and query/subject.length and <tol>/<both.ends> [DEFAULT: TRUE]
6) <tol|NUM>:               alignments ending at <tol> basepairs from a query/subject start/end end will still be considered as being at an alignment end [DEFAULT: 5]
7) <both.ends|BOOLEAN>:      if TRUE, both query and subject alignments need to be at an alignment end ; if FALSE, one is enough [DEFAULT: FALSE]
8) <min.alnlen|NUM>:        minimum alignment length for hit to be considered [DEFAULT: 0]
9) <min.percid|NUM>:        minimum percent identity for hit to be considered [DEFAULT: 0]
```

**Depends on**
```
--
```


**Example**
```
# identify consensus sequence (alignment) overlaps assuming that these two sequences are on the same linkage group (01):
# '>NC_033804.1_LG_01_10292880_10292979_ID_658.merged'
# '>NC_033804.1_LG_01_12382290_12382509_ID_1154'

find.overlapping.alignments.R consensus.vs.self.blast.filtered TRUE 'LG_' '_'
```

## 5) align.overlapping.contigs.sh

**Usage**
```
align.overlapping.contigs.sh -l <file> -c <directory> -m <string> -t <positive integer>
```

**Arguments**
```
# Required
-l                 Path with list of regions to be merged. Usually the output of find.overlapping.alignments.R,
                   which produces a file ending in '.blast.filtered.list'.
-c                 Path to directory with multifasta files, containing FASTA files with unaligned contigs.

# Optional
-m  ['localpair']  Alignment model. Passed on to mafft.
-x  ['']           Prefix of sequence name in input file (-l option) but missing in multifasta sequence names.
-y  ['']           Suffix of sequence name in input file (-l option) but missing in multifasta sequence names.
-t  [4]            Number of parallel threads.
```

**Details**
```
This script creates an output directory of the form 'mafft.overlapping.${N}', where ${N} is the number of 
overlapping alignments passed with -l.
```

**Depends on**
```
mafft
awk
get.overlap.consensus.R
```


**Example**
```
align.overlapping.contigs.sh -l consensus.blast.filtered.list -c multifasta.63.2396 -m 'localpair' -t 20
```

## 6) filter.merged.alignments.sh

**Usage**
```
filter.merged.alignments.sh -d <directory> -s <numeric fraction>
```

**Arguments**
```
# Required
-d          Path to directory with merged alignments.

# Optional
-s  [0.95]  Minimum alignment score [between 0 and 1]. Overlapping alignments will only be retained if the
            merging procedure in `align.overlapping.contigs.sh` was successful, as indicated by the 
            visualizations and/or the alignment score.
```

**Details**
```
This script filters merged alignments by moving unsuccessfully merged alignments to a separate subdirectory.
```

**Depends on**
```
--
```


**Example**
```
filter.merged.alignments.sh -d mafft.overlapping.15 -s 0.95
```


## Apply trimming as in [Step 5](Step5_Alignment_and_Alignment_Trimming.md)


## 7) replace.overlapping.alignments.R

**Usage**
```
replace.overlapping.alignments.R <directory> <directory> <file>
```

**Arguments**
```
# Required
1) <alndir|CHR>:     path to directory with overlapping fasta alignments to be replaced (mafft.NIND.NLOC.Y) ; 
2) <mrgdir|CHR>:     path to directory with corresponding merged (and filtered) fasta alignments ('mafft.overlapping.X.Y') ;
3) <mrglist|CHR>:    path to overlapping loci list ('mafft.NIND.NLOC.Y.consZ.vs.self.blast.filtered.list') ;
       
# Optional
4) <revert|BOOLEAN>: if TRUE, will revert any previous filtering applied [DEFAULT: TRUE]
```

**Details**
```
This script replaces overlapping alignments in the first directory with filtered merged alignments in the
second directory.
```

**Depends on**
```
--
```


**Example**
```
replace.overlapping.alignments.R mafft.63.2396.c0.5.d0.25.c0.4.n0.5 mafft.overlapping.15 consensus.blast.filtered.list
```

## Continue
[➜ Continue to Step 7](Step7_Create_Representative_Reference_Sequences.md)
