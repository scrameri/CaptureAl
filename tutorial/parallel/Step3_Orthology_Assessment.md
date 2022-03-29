[← Back to Step 2](Step2_Sequence_Assembly.md)


# STEP 3

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step3.png)


## 1) select.best.contigs.per.locus.sh

**Usage**
```
select.best.contigs.per.locus.sh -s $s -l $l -r $1 -d $assemblies -t 20
```

**Arguments**
```
# Required
-s  .txt file with names of retained samples (without header or '>')
-l  .txt file with names of retained target regions (without header or '>')
-r  FASTA file with reference sequences
-d  path to input directory with assembled contigs

# Optional
-c  path from input directory to contigs, will replace 'SAMPLE' and 'LOCUS' with the respective strings
    [DEFAULT: `SAMPLE.dipspades/extracted_reads_SAMPLE.fastq.LOCUS.ids.spades/dipspades/consensus_contigs.fasta`]
-g  FLAG, if given, will NOT write query and target ranges to .exonerate output> [DEFAULT: --]
-m  alignment model [DEFAULT: `affine:local`]
-t  number of threads [DEFAULT: 15]


```

**Details**
```
# EXONERATE gapped alignment options
- affine:global
This performs gapped global alignment, similar to the Needleman-Wunsch algorithm, except with affine gaps. Global alignment requires that both the sequences in their entirety are included in the alignment.

- affine:bestfit
This performs a best fit or best location alignment of the query onto the target sequence. The entire query sequence will be included in the alignment, but only the best location for its alignment on the target sequence.

- affine:local
This is local alignment with affine gaps, similar to the Smith-Waterman-Gotoh algorithm. A general-purpose alignment algorithm. As this is local alignment, any subsequence of the query and target sequence may appear in the alignment.
```

**Depends on**
```
exonerate
extract.all.fasta.seqs.R
```


**Example**
```
select.best.contigs.per.locus.sh -s samples.txt -l loci.txt -r reference.fasta -d NovaSeq-run1_assemblies -t 20
```


## 2) combine.contigs.parallel.sh

**Usage**
```
combine.contigs.parallel.sh -s <file> -d <directory> -a <integer> -c <integer> -t <integer>
```

**Arguments**
```
# Required
-s  .txt file with sample names (without header or '>')
-d  path to directory with exonerate results

# Optional
-a  minimum target alignment length [DEFAULT: 80]
-c  minimum normalized alignment score [DEFAULT: 2]
-t  number of threads [DEFAULT: 4]

```

**Depends on**
```
combine.contigs.R
```


**Example**
```
combine.contigs.parallel.sh -s samples.txt -d NovaSeq-run1_exonerate -a 1 -c 1 -t 20
```


## 3) collect.exonerate.stats.R

**Usage**
```
collect.exonerate.stats.R <sample file> <OPTIONAL: directory> <OPTIONAL: string>
```

**Arguments**
```
# Required
1) <sfile|CHR>:  path to sample file containing paths to mapping directories

# Optional (if one is given, all must be given in this order)
2) <dir|CHR>:    folder with subdirectories for each sample [DEFAULT: current directory ]
3) <suffix|CHR>: sample suffix present in <sfile> and absent in <dir> subdirectories  [DEFAULT: '']
```

**Example**
```
collect.exonerate.stats.R samples.txt NovaSeq-run1_exonerate
```


## 4) plot.contig.numbers.R

**Usage**
```
plot.contig.numbers.R <contig number file> <taxon group file>
```

**Arguments**
```
# Required
1) <file|CHR>: path to collected contig numbers (loci_contignumbers.txt)

# Optional
2) <meta|CHR>: path to metadata file mapping individuals (1st column) to groups (2nd column). Header and tab separation expected. More individuals in different order than in <file> are ok.
```

**Example**
```
plot.contig.numbers.R loci_contignumbers.txt samples.mapfile.txt
```


## Continue
[➜ Continue to Step 4](Step4_Sample_and_Locus_Filtering.md)
