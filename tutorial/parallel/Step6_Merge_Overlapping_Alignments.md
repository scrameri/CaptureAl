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


## 3) [blast.vs.self.sh](https://github.com/scrameri/CaptureAl/wiki/blast.vs.self.sh)

This determines [BLAST+](https://scicomp.ethz.ch/public/manual/BLAST/BLAST.pdf) hits between pairs of different consensus sequences (i.e., alignments), and writes the filtered results (keeping one out of two reciprocal hits) to a file `*.blast.vs.self.blast.filtered`.

**Example**
```
blast.vs.self.sh consensus.fasta
```

## 4) [find.overlapping.alignments.R](https://github.com/scrameri/CaptureAl/wiki/find.overlapping.alignments.R)

This identifies target regions with physical overlap based on [BLAST+](https://scicomp.ethz.ch/public/manual/BLAST/BLAST.pdf) hits among alignment consensus sequences.

**Example**
```
# identify consensus sequence (alignment) overlaps assuming that these two sequences are on the same linkage group (01):
# '>NC_033804.1_LG_01_10292880_10292979_ID_658.merged'
# '>NC_033804.1_LG_01_12382290_12382509_ID_1154'

find.overlapping.alignments.R consensus.vs.self.blast.filtered TRUE 'LG_' '_'
```

## 5) [align.overlapping.contigs.sh](https://github.com/scrameri/CaptureAl/wiki/align.overlapping.contigs.sh)

This performs an alignment of contigs from overlapping target regions listed in `consensus.blast.filtered.list`.

An output directory is created named 'mafft.overlapping.${N}', where ${N} is the number of overlapping alignments passed with `-l`.

That directory contains the merged alignment (to be used for trimming and to be combined with alignments of non-overlapping regions), but also statistics (`*.log`) and a visualization (`*.pdf`) of the merging procedure, which can be used to decide whether merging was successful.

**Example**
```
align.overlapping.contigs.sh -l consensus.blast.filtered.list -c multifasta.63.2396 -m 'localpair' -t 20
```

## 6) [filter.merged.alignments.sh](https://github.com/scrameri/CaptureAl/wiki/filter.merged.alignments.sh)

This filters any overlapping alignments in `mafft.overlapping.15` that have an alignment score < 0.95.

To decide on the optimal value for `-s`, it is recommended to compare the alignment scores (`*.log`) with the alignment visualiztions (`*.pdf`).

**Example**
```
filter.merged.alignments.sh -d mafft.overlapping.15 -s 0.95
```


## Apply trimming as in [Step 5](Step5_Alignment_and_Alignment_Trimming.md)


## 7) [replace.overlapping.alignments.R](https://github.com/scrameri/CaptureAl/wiki/replace.overlapping.alignments.R)

This replaces two or more (trimmed) overlapping alignments in `mafft.63.2396.c0.5.d0.25.c0.4.n0.5` with a corresponding single merged (trimmed) alignment, for all successfully merged alignments.

**Example**
```
replace.overlapping.alignments.R mafft.63.2396.c0.5.d0.25.c0.4.n0.5 mafft.overlapping.15 consensus.blast.filtered.list
```

## Continue
[➜ Continue to Step 7](Step7_Create_Representative_Reference_Sequences.md)
