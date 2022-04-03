[← Back to Step 4](Step4_Sample_and_Locus_Filtering.md) ................................................................................................................................................................. [➜ Continue to Step 6](Step6_Merge_Overlapping_Alignments.md)


# STEP 5

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step5.png)


## 1) [create.multifastas.parallel.sh](https://github.com/scrameri/CaptureAl/wiki/create.multifastas.parallel.sh)

This creates a multi-sequence FASTA file with unaligned contigs of multiple samples specified in [samples.txt](https://github.com/scrameri/CaptureAl/blob/master/tutorial/data/samples.txt), using contigs located in sample subdirectories of `NovaSeq-run1_exonerate`, for a batch of target regions specified in [loci.txt](https://github.com/scrameri/CaptureAl/blob/master/tutorial/data/loci.txt) in parallel.

**Example**
```
create.multifastas.parallel.sh -s samlples.txt -l loci.txt -d NovaSeq-run1_exonerate -t 20
```


## 2) [align.multifastas.parallel.sh](https://github.com/scrameri/CaptureAl/wiki/align.multifastas.parallel.sh)

This uses [mafft](https://mafft.cbrc.jp/alignment/software/manual/manual.html) to alignm multifasta sequences located in the `multifasta.63.2396` directory in parallel. Consult the [mafft](https://mafft.cbrc.jp/alignment/software/manual/manual.html) manpage for the choice of alignment model (`-m` option).

**Example**
```
align.multifastas.parallel.sh -d multifasta.63.2396 -m 'localpair' -t 20
```


## 3) [trim.alignment.ends.parallel.sh](https://github.com/scrameri/CaptureAl/wiki/trim.alignment.ends.parallel.sh)

This uses the raw alignments in the `mafft.63.2396` directory and creates a new output directory with trimmed alignments.

The [mapfile.txt](https://github.com/scrameri/CaptureAl/blob/master/tutorial/data/mapfile.txt) is expected to have a header, and the first column is interpreted to contain the identifiers of samples to be included in trimmed alignments.

**Example**
```
# no visualization
trim.alignment.ends.parallel.sh -s mapfile.txt -d mafft.63.2396 -c 0.5 -n 0.25 -t 20

# with visualization
trim.alignment.ends.parallel.sh -s mapfile.txt -d mafft.63.2396 -c 0.5 -n 0.25 -t 20 -v
```

## 4) [trim.alignments.parallel.sh](https://github.com/scrameri/CaptureAl/wiki/trim.alignments.parallel.sh)

This uses the end-trimmed alignments in the `mafft.63.2396.c0.5.n0.25` directory and creates a new output directory with internally trimmed alignments.

The [mapfile.txt](https://github.com/scrameri/CaptureAl/blob/master/tutorial/data/mapfile.txt) is expected to have a header, and the first column is interpreted to contain the identifiers of samples to be included in trimmed alignments.

**Example**
```
# no visualization
trim.alignments.parallel.sh -s mapfile.txt -d mafft.63.2396.c0.5.n0.25 -c 0.4 -z 20 -S 1 -n 0.5 -t 20

# with visualization
trim.alignments.parallel.sh -s mapfile.txt -d mafft.63.2396.c0.5.n0.25 -c 0.4 -z 20 -S 1 -n 0.5 -t 20 -v

# with visualization and internal trimming
trim.alignments.parallel.sh -s mapfile.txt -d mafft.63.2396.c0.5.n0.25 -c 0.4 -z 20 -S 1 -n 0.5 -t 20 -iv
```

## Continue
[← Back to Step 4](Step4_Sample_and_Locus_Filtering.md) ................................................................................................................................................................. [➜ Continue to Step 6](Step6_Merge_Overlapping_Alignments.md)
